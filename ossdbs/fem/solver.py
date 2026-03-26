# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from abc import ABC, abstractmethod

import ngsolve
import numpy as np

from ossdbs.fem.preconditioner import (
    BDDCPreconditioner,
    JacobiPreconditioner,
    Preconditioner,
)

_logger = logging.getLogger(__name__)


def _is_finite_scalar(value) -> bool:
    """Return whether a scalar or complex scalar is finite."""
    return np.isfinite(value.real) and np.isfinite(value.imag)


def _test_matrix_on_freedofs(
    matrix, template_vector, freedofs, ntests: int = 3
) -> None:
    """Log symmetry and SPD checks on the free DOFs only."""
    projector = ngsolve.Projector(freedofs, True)

    _logger.debug("Testing matrix symmetry/SPD on free DOFs")
    for i in range(ntests):
        x = template_vector.CreateVector()
        y = template_vector.CreateVector()
        x.SetRandom()
        y.SetRandom()
        x.data = projector * x
        y.data = projector * y

        ax = template_vector.CreateVector()
        ay = template_vector.CreateVector()
        ax.data = projector * (matrix * x)
        ay.data = projector * (matrix * y)

        spd_val = x.InnerProduct(ax)
        sym_lhs = x.InnerProduct(ay)
        sym_rhs = ax.InnerProduct(y)
        sym_diff = abs(sym_lhs - sym_rhs)

        _logger.debug(
            "free-dof test %s: <x, A x> = %s, symmetry diff = %s",
            i,
            spd_val,
            sym_diff,
        )

        if not _is_finite_scalar(spd_val) or not _is_finite_scalar(sym_diff):
            _logger.debug("free-dof test %s: matrix test produced nan/inf", i)
            return
        if spd_val.real <= 0:
            _logger.debug("free-dof test %s: matrix is not positive definite", i)
            return


def _build_jacobi_inverse_diag(matrix, freedofs) -> np.ndarray:
    """Build inverse diagonal for a Jacobi preconditioner on free DOFs only."""
    rows, cols, vals = matrix.COO()
    rows = np.array(rows, dtype=int)
    cols = np.array(cols, dtype=int)
    vals = np.array(vals)
    freedofs_np = np.array(freedofs, dtype=bool)

    diag = np.zeros(matrix.height, dtype=vals.dtype)
    diagonal_entries = rows == cols
    diag[rows[diagonal_entries]] = vals[diagonal_entries]

    inv_diag = np.zeros_like(diag)
    valid = freedofs_np & (~np.isclose(diag, 0.0))
    inv_diag[valid] = 1.0 / diag[valid]
    return inv_diag


def _apply_jacobi_preconditioner(dst, src, inv_diag: np.ndarray, freedofs) -> None:
    """Apply a diagonal Jacobi inverse to a vector on the free DOFs only."""
    src_np = src.FV().NumPy()
    dst_np = dst.FV().NumPy()
    dst_np[:] = inv_diag * src_np
    dst_np[~np.array(freedofs, dtype=bool)] = 0.0


def _create_customized_local_preconditioner(
    matrix,
    freedofs,
    template_vec,
    residual=None,
    debug_enabled: bool = False,
) -> JacobiPreconditioner:
    """Create a custom Jacobi preconditioner and optionally log diagnostics."""
    inv_diag = _build_jacobi_inverse_diag(matrix, freedofs)
    if debug_enabled:
        _logger.debug(
            "Jacobi inverse diagonal: nonzero=%s finite=%s",
            np.count_nonzero(inv_diag),
            np.count_nonzero(np.isfinite(inv_diag)),
        )
        if residual is not None:
            z = residual.CreateVector()
            _apply_jacobi_preconditioner(z, residual, inv_diag, freedofs)
            _logger.debug("Norm(M^-1*r) = %s", z.Norm())
            _logger.debug("<M^-1*r, r> = %s", z.InnerProduct(residual))
    return JacobiPreconditioner(inv_diag, freedofs, template_vec)


def _get_debug_mode() -> tuple[bool, bool]:
    """Return whether debug logging is enabled and whether rates should print."""
    debug_enabled = _logger.isEnabledFor(logging.DEBUG)
    printrates = debug_enabled
    return debug_enabled, printrates


def _log_initial_system_state(
    bilinear_form: ngsolve.BilinearForm,
    linear_form: ngsolve.LinearForm,
    grid_function: ngsolve.GridFunction,
    debug_enabled: bool,
    solver_name: str,
):
    """Build and optionally log the initial residual information."""
    if debug_enabled:
        _logger.debug("Linear form norm = %s", linear_form.vec.Norm())
        _test_matrix_on_freedofs(
            bilinear_form.mat,
            grid_function.vec,
            grid_function.space.FreeDofs(),
        )

    residual = linear_form.vec.CreateVector()
    residual.data = linear_form.vec - bilinear_form.mat * grid_function.vec
    if debug_enabled:
        _logger.debug("Residual norm before %s = %s", solver_name, residual.Norm())
        tmp_a = linear_form.vec.CreateVector()
        tmp_a.data = bilinear_form.mat * residual
        _logger.debug("Norm(A*r) = %s", tmp_a.Norm())
        _logger.debug("<r, A*r> = %s", residual.InnerProduct(tmp_a))
    return residual


def _log_native_preconditioner_state(
    preconditioner, residual, debug_enabled: bool
) -> None:
    """Log the native preconditioner action on the initial residual."""
    if not debug_enabled:
        return
    tmp_pre = residual.CreateVector()
    tmp_pre.data = preconditioner * residual
    _logger.debug("Norm(pre*r) = %s", tmp_pre.Norm())
    _logger.debug("<pre*r, r> = %s", tmp_pre.InnerProduct(residual))


def _warn_local_preconditioner_issue(
    precond_type: str,
    iterations: int,
    residual,
    correction,
    final_residual,
    debug_enabled: bool,
) -> None:
    """Warn when the native local preconditioner exits suspiciously fast."""
    if precond_type != "local":
        return

    suspicious_fast_exit = iterations <= 3
    suspicious_zero_correction = residual.Norm() > 0 and correction.Norm() == 0
    suspicious_final_residual = False
    if debug_enabled and final_residual is not None:
        suspicious_final_residual = (
            not np.isfinite(final_residual.Norm())
            or final_residual.Norm() >= residual.Norm()
        )
    if suspicious_fast_exit and (
        suspicious_zero_correction or suspicious_final_residual
    ):
        _logger.warning(
            "Native NGSolve 'local' preconditioner may be unstable for this "
            "problem. If the result looks suspicious, try changing "
            '"Preconditioner": "local" to "customized_local" in the JSON input.'
        )


def _finalize_krylov_solve(
    solver,
    bilinear_form: ngsolve.BilinearForm,
    linear_form: ngsolve.LinearForm,
    grid_function: ngsolve.GridFunction,
    residual,
    correction,
    debug_enabled: bool,
    solver_name: str,
    maxsteps: int,
    precond_type: str | None = None,
) -> None:
    """Apply the correction, log final diagnostics, and emit warnings."""
    grid_function.vec.data += correction

    residual_after = None
    if debug_enabled:
        residual_after = linear_form.vec.CreateVector()
        residual_after.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        _logger.debug("Converged after %s iterations", solver.iterations)
        if solver.residuals:
            _logger.debug("Final %s residual = %s", solver_name, solver.residuals[-1])
        _logger.debug("Residual norm after %s = %s", solver_name, residual_after.Norm())

    _warn_local_preconditioner_issue(
        precond_type,
        solver.iterations,
        residual,
        correction,
        residual_after,
        debug_enabled,
    )

    if solver.iterations >= maxsteps:
        _logger.warning(
            "Did not converge after %s iterations! "
            "Increase the maximum number of steps!",
            solver.iterations,
        )


class Solver(ABC):
    """Computes the boundary volume problem."""

    DEFAULT_PRECONDITIONER = BDDCPreconditioner()

    def __init__(
        self,
        precond_par: Preconditioner = DEFAULT_PRECONDITIONER,
        maxsteps: int = 10000,
        precision: float = 1e-12,
    ) -> None:
        """Initialize the solver.

        Parameters
        ----------
        precond_par : Preconditioner
            Preconditioner
        maxsteps : int
            Maximum steps before solver ends
        precision : float
            Desired precision
        """
        self._precond_par = precond_par.to_dictionary()
        self._maxsteps = maxsteps
        self._precision = precision

    @abstractmethod
    def bvp(
        self,
        bilinear_form: ngsolve.BilinearForm,
        linear_form: ngsolve.LinearForm,
        grid_function: ngsolve.GridFunction,
    ) -> None:
        """Solve boundary-value problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
            Bilinear form (left-hand side)
        linear_form: ngsolve.LinearForm
            Linear form (right-hand side)
        grid_function: ngsolve.GridFunction
            Solution vector

        """
        pass


class CGSolver(Solver):
    """Conjugate gradient solver."""

    def bvp(
        self,
        bilinear_form: ngsolve.BilinearForm,
        linear_form: ngsolve.LinearForm,
        grid_function: ngsolve.GridFunction,
    ) -> None:
        """Solve boundary-value problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
            Bilinear form (left-hand side)
        linear_form: ngsolve.LinearForm
            Linear form (right-hand side)
        grid_function: ngsolve.GridFunction
            Solution vector

        """
        preconditioner = None
        if self._precond_par.get("type") != "customized_local":
            preconditioner = ngsolve.Preconditioner(
                bf=bilinear_form, **self._precond_par
            )
        bilinear_form.Assemble()
        linear_form.Assemble()
        debug_enabled, printrates = _get_debug_mode()
        residual = _log_initial_system_state(
            bilinear_form,
            linear_form,
            grid_function,
            debug_enabled,
            "CG",
        )

        if self._precond_par.get("type") == "customized_local":
            preconditioner = _create_customized_local_preconditioner(
                bilinear_form.mat,
                grid_function.space.FreeDofs(),
                grid_function.vec,
                residual=residual,
                debug_enabled=debug_enabled,
            )
        else:
            _log_native_preconditioner_state(preconditioner, residual, debug_enabled)
        solver = ngsolve.krylovspace.CGSolver(
            bilinear_form.mat,
            pre=preconditioner,
            printrates=printrates,
            maxiter=self._maxsteps,
            tol=self._precision,
        )
        corr = grid_function.vec.CreateVector()
        corr[:] = 0.0
        if debug_enabled:
            correction_label = (
                "Jacobi-CG"
                if self._precond_par.get("type") == "customized_local"
                else "CG"
            )
            solver.Solve(rhs=residual, sol=corr, initialize=True)
            _logger.debug(
                "Correction norm after %s = %s", correction_label, corr.Norm()
            )
        else:
            solver.Solve(rhs=residual, sol=corr, initialize=True)
        _finalize_krylov_solve(
            solver,
            bilinear_form,
            linear_form,
            grid_function,
            residual,
            corr,
            debug_enabled,
            "CG",
            self._maxsteps,
            self._precond_par.get("type"),
        )


class GMRESSolver(Solver):
    """GMRes solver."""

    def bvp(
        self,
        bilinear_form: ngsolve.BilinearForm,
        linear_form: ngsolve.LinearForm,
        grid_function: ngsolve.GridFunction,
    ) -> None:
        """Solve boundary-value problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
            Bilinear form (left-hand side)
        linear_form: ngsolve.LinearForm
            Linear form (right-hand side)
        grid_function: ngsolve.GridFunction
            Solution vector

        """
        preconditioner = None
        if self._precond_par.get("type") != "customized_local":
            preconditioner = ngsolve.Preconditioner(
                bf=bilinear_form, **self._precond_par
            )
        bilinear_form.Assemble()
        linear_form.Assemble()
        debug_enabled, printrates = _get_debug_mode()
        residual = _log_initial_system_state(
            bilinear_form,
            linear_form,
            grid_function,
            debug_enabled,
            "GMRES",
        )
        if self._precond_par.get("type") == "customized_local":
            preconditioner = _create_customized_local_preconditioner(
                bilinear_form.mat,
                grid_function.space.FreeDofs(),
                grid_function.vec,
                residual=residual,
                debug_enabled=debug_enabled,
            )
        else:
            _log_native_preconditioner_state(preconditioner, residual, debug_enabled)
        solver = ngsolve.krylovspace.GMResSolver(
            mat=bilinear_form.mat,
            pre=preconditioner,
            printrates=printrates,
            maxiter=self._maxsteps,
            tol=self._precision,
        )
        corr = grid_function.vec.CreateVector()
        corr[:] = 0.0
        if debug_enabled:
            solver.Solve(rhs=residual, sol=corr, initialize=True)
            _logger.debug("Correction norm after GMRES = %s", corr.Norm())
        else:
            solver.Solve(rhs=residual, sol=corr, initialize=True)
        _finalize_krylov_solve(
            solver,
            bilinear_form,
            linear_form,
            grid_function,
            residual,
            corr,
            debug_enabled,
            "GMRES",
            self._maxsteps,
            self._precond_par.get("type"),
        )


class DirectSolver(Solver):
    """Direct solver."""

    def bvp(
        self,
        bilinear_form: ngsolve.BilinearForm,
        linear_form: ngsolve.LinearForm,
        grid_function: ngsolve.GridFunction,
    ) -> None:
        """Solve boundary-value problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
            Bilinear form (left-hand side)
        linear_form: ngsolve.LinearForm
            Linear form (right-hand side)
        grid_function: ngsolve.GridFunction
            Solution vector

        Notes
        -----
        TODO have a second look
        """
        bilinear_form.Assemble()
        linear_form.Assemble()
        inverse = bilinear_form.mat.Inverse(
            freedofs=grid_function.space.FreeDofs(), inverse="pardiso"
        )
        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r
