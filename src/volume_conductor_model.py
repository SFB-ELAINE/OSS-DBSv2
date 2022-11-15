
from src.mesh import Mesh
from src.laplace_equation import LaplaceEquation
import ngsolve

from src.voxel_space import VoxelSpace


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
                 conductivity: VoxelSpace,
                 permitivity: VoxelSpace = None,
                 ) -> None:

        self.__conductivity = conductivity
        self.__permitivity = permitivity
        self.__complex = True if self.__permitivity else False

    def evaluate_potential(self,
                           mesh: Mesh = None,
                           boundaries: dict = None) \
            -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Returns
        -------

        return : tuple
            Postprocessed data: lectric_field, V_contact, Power, potential

        """
        start, end = mesh.bounding_box()
        type = 'complex' if self.__permitivity else 'float'
        conductivity = self.__conductivity.data.astype(type)
        sigma = ngsolve.VoxelCoefficient(start=start,
                                         end=end,
                                         values=conductivity,
                                         linear=False)

        if self.__permitivity:
            permitivity = self.__permitivity.data.astype(type)
            permitivities = ngsolve.VoxelCoefficient(start=start,
                                                     end=end,
                                                     values=permitivity,
                                                     linear=False)
            sigma += permitivities

        space = mesh.sobolev_space(complex=self.__complex)
        potential = ngsolve.GridFunction(space=space)
        coefficient = mesh.boundary_coefficients(boundaries=boundaries)

        potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=space, coefficient=sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential, self.__error(potential, mesh)

    def __error(self, potential: ngsolve.comp.GridFunction, mesh: Mesh) \
            -> ngsolve.fem.CoefficientFunction:
        flux = ngsolve.grad(potential)
        space = mesh.flux_space(complex=self.__complex)
        flux_potential = ngsolve.GridFunction(space=space)
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)
