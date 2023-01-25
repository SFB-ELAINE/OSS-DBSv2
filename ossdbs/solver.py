import ngsolve

from ossdbs.preconditioner import BDDCPreconditioner, PreconditionerParameters
from abc import ABC, abstractmethod


class Solver(ABC):

    @abstractmethod
    def bvp(self, bilinear_form, linear_form, grid_function):
        pass


class CGSolver:

    def __init__(self,
                 precond_par: PreconditionerParameters = BDDCPreconditioner(),
                 printrates: bool = True,
                 maxsteps: int = 10000,
                 precision: float = 1e-12) -> None:
        self.__precond_par = precond_par.to_dictionary()
        self.__printrates = printrates
        self.__maxsteps = maxsteps
        self.__precision = precision

    def bvp(self, bilinear_form, linear_form, grid_function):
        preconditioner = ngsolve.Preconditioner(bf=bilinear_form,
                                                **self.__precond_par)

        bilinear_form.Assemble()
        linear_form.Assemble()
        inverse = ngsolve.CGSolver(mat=bilinear_form.mat,
                                   pre=preconditioner.mat,
                                   printrates=self.__printrates,
                                   maxsteps=self.__maxsteps,
                                   precision=self.__precision)

        r = linear_form.vec.CreateVector()
        r.data = linear_form.vec - bilinear_form.mat * grid_function.vec
        grid_function.vec.data = grid_function.vec.data + inverse * r


class GMRESSolver:

    def __init__(self,
                 precond_par: PreconditionerParameters = BDDCPreconditioner(),
                 printrates: bool = True,
                 maxsteps: int = 10000,
                 precision: float = 1e-12) -> None:
        self.__precond_par = precond_par.to_dictionary()
        self.__printrates = printrates
        self.__maxsteps = maxsteps
        self.__precision = precision

    def bvp(self, bilinear_form, linear_form, grid_function):
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
