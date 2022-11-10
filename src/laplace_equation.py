import ngsolve
import numpy as np

class LaplaceEquation:

    def __init__(self,
                 space: ngsolve.comp.H1,
                 coefficient: ngsolve.fem.CoefficientFunction) -> None:

        u = space.TrialFunction()
        v = space.TestFunction()
        self.__a = ngsolve.BilinearForm(space=space, symmetric=True)
        self.__a += (coefficient *
                     ngsolve.grad(u) *
                     ngsolve.grad(v) *
                     ngsolve.dx)
        self.__f = ngsolve.LinearForm(space=space)
        self.__preconditioner = ngsolve.Preconditioner(bf=self.__a,
                                                       type="bddc",
                                                       coarsetype="h1amg")

    def solve_bvp(self, input: ngsolve.comp.GridFunction) \
            -> ngsolve.la.DynamicVectorExpression:
        """Solve boundary value problem."""
        self.__a.Assemble()
        self.__f.Assemble()
        inverse = ngsolve.CGSolver(mat=self.__a.mat,
                                   pre=self.__preconditioner.mat,
                                   printrates=True,
                                   maxsteps=10000,
                                   precision=1e-12)
        r = self.__f.vec.CreateVector()
        r.data = self.__f.vec - self.__a.mat * input.vec
        return input.vec.data + inverse * r
