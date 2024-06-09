# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from abc import ABC, abstractmethod

import ngsolve

from ossdbs.fem.preconditioner import BDDCPreconditioner, Preconditioner

_logger = logging.getLogger(__name__)


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
        preconditioner = ngsolve.Preconditioner(bf=bilinear_form, **self._precond_par)

        bilinear_form.Assemble()
        linear_form.Assemble()
        inverse = ngsolve.CGSolver(
            mat=bilinear_form.mat,
            pre=preconditioner.mat,
            printrates=False,
            maxsteps=self._maxsteps,
            precision=self._precision,
        )
        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r

        # check the number of iterations
        _logger.debug(f"Converged after {inverse.GetSteps()} iterations")
        if inverse.GetSteps() >= self._maxsteps:
            _logger.warning(
                f"Did not converge after {inverse.GetSteps()} iterations!"
                " Increase the maximum number of steps!"
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
        preconditioner = ngsolve.Preconditioner(bf=bilinear_form, **self._precond_par)

        bilinear_form.Assemble()
        linear_form.Assemble()
        inverse = ngsolve.GMRESSolver(
            mat=bilinear_form.mat,
            pre=preconditioner.mat,
            printrates=False,
            maxsteps=self._maxsteps,
            precision=self._precision,
        )

        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r

        # check the number of iterations
        _logger.debug(f"Converged after {inverse.GetSteps()} iterations")
        # tested it and apparently the steps are slightly miscounted
        if inverse.GetSteps() >= self._maxsteps - 2:
            _logger.warning(
                f"Did not converge after {inverse.GetSteps()} iterations!"
                " Increase the maximum number of steps!"
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
