
from src.mesh import Mesh
from src.laplace_equation import LaplaceEquation
import ngsolve


class VolumeConductor:
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

    def __init__(self,
                 conductivity: ngsolve.CoefficientFunction,
                 ) -> None:

        self.__sigma = conductivity

    def evaluate_potential(self, mesh: Mesh = None) \
            -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Returns
        -------

        return : tuple
            Postprocessed data: lectric_field, V_contact, Power, potential

        """

        space = mesh.sobolev_space()
        potential = ngsolve.GridFunction(space=space)
        coefficient = mesh.boundary_coefficients()
        potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=space, coefficient=self.__sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential, self.__error(potential, mesh)

    def __error(self, potential: ngsolve.comp.GridFunction, mesh: Mesh) \
            -> ngsolve.fem.CoefficientFunction:
        flux = ngsolve.grad(potential)
        flux_potential = ngsolve.GridFunction(space=mesh.flux_space())
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)
