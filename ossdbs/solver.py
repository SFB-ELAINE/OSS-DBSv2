import ngsolve

from ossdbs.preconditioner import BDDCPreconditioner, Preconditioner
from abc import ABC, abstractmethod


class Solver(ABC):
    """Computes the boundary volume problem.
    """

    @abstractmethod
    def bvp(self,
            bilinear_form: ngsolve.BilinearForm,
            linear_form: ngsolve.LinearForm,
            grid_function: ngsolve.GridFunction) -> None:
        pass


class CGSolver(Solver):
    """

    Parameters
    ----------
    precond_par : Preconditioner
    printrates : bool
    maxsteps : int
    precision : float
    """
    def __init__(self,
                 precond_par: Preconditioner = BDDCPreconditioner(),
                 printrates: bool = True,
                 maxsteps: int = 10000,
                 precision: float = 1e-12) -> None:
        self.__precond_par = precond_par.to_dictionary()
        self.__printrates = printrates
        self.__maxsteps = maxsteps
        self.__precision = precision

    def bvp(self,
            bilinear_form: ngsolve.BilinearForm,
            linear_form: ngsolve.LinearForm,
            grid_function: ngsolve.GridFunction) -> None:
        """Compute boundary volume problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
        linear_form: ngsolve.LinearForm
        grid_function: ngsolve.GridFunction
        """
        print("preconditioner")
        preconditioner = ngsolve.Preconditioner(bf=bilinear_form,
                                                **self.__precond_par)
        print("assemble:")
        bilinear_form.Assemble()
        linear_form.Assemble()
        print("solve:")
        inverse = ngsolve.CGSolver(mat=bilinear_form.mat,
                                   pre=preconditioner.mat,
                                   printrates=self.__printrates,
                                   maxsteps=self.__maxsteps,
                                   precision=self.__precision)

        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r


class GMRESSolver(Solver):
    """

    Parameters
    ----------
    precond_par : Preconditioner
    printrates : bool
    maxsteps : int
    precision : float
    """
    def __init__(self,
                 precond_par: Preconditioner = BDDCPreconditioner(),
                 printrates: bool = True,
                 maxsteps: int = 10000,
                 precision: float = 1e-12) -> None:
        self.__precond_par = precond_par.to_dictionary()
        self.__printrates = printrates
        self.__maxsteps = maxsteps
        self.__precision = precision

    def bvp(self,
            bilinear_form: ngsolve.BilinearForm,
            linear_form: ngsolve.LinearForm,
            grid_function: ngsolve.GridFunction) -> None:
        """Compute boundary volume problem.

        Parameters
        ----------
        bilinear_form: ngsolve.BilinearForm
        linear_form: ngsolve.LinearForm
        grid_function: ngsolve.GridFunction
        """
        preconditioner = ngsolve.Preconditioner(bf=bilinear_form,
                                                **self.__precond_par)

        bilinear_form.Assemble()
        linear_form.Assemble()
        inverse = ngsolve.GMRESSolver(mat=bilinear_form.mat,
                                      pre=preconditioner.mat,
                                      printrates=self.__printrates,
                                      maxsteps=self.__maxsteps,
                                      precision=self.__precision)

        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r
