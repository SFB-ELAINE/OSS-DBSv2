
from abc import ABC, abstractmethod
from src.mesh import Mesh
from src.voxels import Voxels
import ngsolve
import numpy as np


class Potential(ABC):

    def __init__(self,
                 gridfunction: ngsolve.GridFunction,
                 mesh: ngsolve.Mesh,
                 ) -> None:
        self._gridfunction = gridfunction
        self._mesh = mesh

    @abstractmethod
    def error(self):
        pass

    def __add__(self, other) -> 'Potential':
        self._gridfunction.vec.data += other._gridfunction.vec.data

    def __mul__(self, value) -> 'Potential':
        self._gridfunction.vec.data *= value


class PotentialEQ:

    def error(self):
        flux = ngsolve.grad(self.__gridfunction)
        space = self.__mesh.flux_space()
        flux_potential = ngsolve.GridFunction(space=space)
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)


class PotentialEQS:

    def error(self):
        flux = ngsolve.grad(self.__gridfunction)
        space = self.__mesh.flux_space(complex=True)
        flux_potential = ngsolve.GridFunction(space=space)
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)


class VolumeConductorQS:
    """Model for representing a volume conductor.

    Attributes
    ----------
    mesh : Mesh

    conductivity : dict

    Methods
    -------
    evaluate_potential(mesh: Mesh, conductivity: dict)
        Evaluate the electric potential of volume conductor.

    """

    def __init__(self, mesh: Mesh, conductivity: Voxels) -> None:
        self.__conductivity = conductivity
        self.__mesh = mesh

    def evaluate_potential(self, boundaries: dict) \
            -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Returns
        -------

        return : tuple
            potential, error

        """

        values = np.real(self.__conductivity.data)
        start, end = self.__conductivity.start, self.__conductivity.end
        sigma = ngsolve.VoxelCoefficient(start, end, values, linear=False)
        space = self.__mesh.sobolev_space()
        potential = ngsolve.GridFunction(space=space)
        coefficient = self.__mesh.boundary_coefficients(boundaries=boundaries)

        potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=space, coefficient=sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential


class VolumeConductorEQS:
    """Model for representing a volume conductor.

    Attributes
    ----------
    mesh : Mesh

    conductivity : dict

    Methods
    -------
    evaluate_potential(mesh: Mesh, conductivity: dict)
        Evaluate the electric potential of volume conductor.

    """

    def __init__(self, mesh: Mesh, conductivity: Voxels) -> None:
        self.__conductivity = conductivity
        self.__mesh = mesh

    def evaluate_potential(self, boundaries: dict) \
            -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Returns
        -------

        return : tuple
            potential, error

        """
        start, end = self.__conductivity.start, self.__conductivity.end
        values = self.__conductivity.data
        sigma = ngsolve.VoxelCoefficient(start, end, values, linear=False)
        space = self.__mesh.sobolev_space(complex=True)
        potential = ngsolve.GridFunction(space=space)
        coefficient = self.__mesh.boundary_coefficients(boundaries=boundaries)

        potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=space, coefficient=sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential


class LaplaceEquation:

    def __init__(self,
                 space: ngsolve.comp.H1,
                 coefficient: ngsolve.fem.CoefficientFunction) -> None:

        u = space.TrialFunction()
        v = space.TestFunction()
        equation = coefficient * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        self.__a = ngsolve.BilinearForm(space=space, symmetric=True)
        self.__a += equation
        self.__f = ngsolve.LinearForm(space=space)
        self.__preconditioner = ngsolve.Preconditioner(bf=self.__a,
                                                       type="bddc",
                                                       coarsetype="local")

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
