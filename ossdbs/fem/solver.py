from abc import ABC, abstractmethod

import ngsolve

from ossdbs.fem.preconditioner import BDDCPreconditioner, Preconditioner


class Solver(ABC):
    """Computes the boundary volume problem."""

    def __init__(
        self,
        precond_par: Preconditioner = BDDCPreconditioner(),
        printrates: bool = True,
        maxsteps: int = 10000,
        precision: float = 1e-12,
    ) -> None:
        """Initialize the solver.

        Parameters
        ----------
        precond_par : Preconditioner
            Preconditioner
        printrates : bool
            Verbose output
        maxsteps : int
            Maximum steps before solver ends
        precision : float
            Desired precision
        """
        self._precond_par = precond_par.to_dictionary()
        self._printrates = printrates
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
            printrates=self._printrates,
            maxsteps=self._maxsteps,
            precision=self._precision,
        )
        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r


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
            printrates=self._printrates,
            maxsteps=self._maxsteps,
            precision=self._precision,
        )

        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r


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
