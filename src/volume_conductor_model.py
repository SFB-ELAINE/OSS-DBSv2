
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

    def evaluate_potential(self, mesh: Mesh = None) -> tuple:
        """Evaluate electrical potential of volume conductor.

        Returns
        -------

        return : tuple
            Postprocessed data: lectric_field, V_contact, Power, potential

        """

        potential = self.__solve_bvp(mesh=mesh)
        return self.__postprocess(potential=potential, mesh=mesh)

    def evaluate_potential2(self,
                            n_dof_limit: int = 1e6,
                            n_refinements: int = 1,
                            mesh: Mesh = None) -> tuple:
        """Evaluate electrical potential of volume conductor.

        Parameters
        ----------
        n_dof_limit : int
            Maximal number of Demensions of Freedom of FEM space.

        n_refinements : int
            Number of mesh refinements.

        Returns
        -------

        return : tuple
            Postprocessed data: lectric_field, V_contact, Power, potential

        """
        potential = self.__solve_bvp(mesh=mesh)
        iterations = 0
        while (mesh.sobolev_space().ndof < n_dof_limit and
               iterations < n_refinements):
            mesh.mark_elements_by_error(error=self.error(potential, mesh))
            mesh.refine()
            potential = self.__solve_bvp(mesh=mesh)
            iterations = iterations + 1
        return self.__postprocess(potential=potential, mesh=mesh)

    def __solve_bvp(self, mesh: Mesh = None) -> ngsolve.comp.GridFunction:
        space = mesh.sobolev_space()
        potential = ngsolve.GridFunction(space=space)
        potential.Set(coefficient=mesh.boundary_coefficients(),
                      VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=space,
                                   coefficient=self.__sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential

    def error(self, potential: ngsolve.comp.GridFunction, mesh: Mesh) \
            -> ngsolve.fem.CoefficientFunction:
        flux = ngsolve.grad(potential)
        flux_potential = ngsolve.GridFunction(space=mesh.flux_space())
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)

    def __postprocess(self, potential: ngsolve.comp.GridFunction, mesh) -> tuple:
        # MappedIntegrationPoint
        electric_field = ngsolve.sqrt(ngsolve.grad(potential) *
                                      ngsolve.grad(potential))
        electric_field = electric_field(
                            mesh.ngsolvemesh()(0, 0, 0.02 / 2)) * 1e3

        V_contact = ngsolve.Integrate(potential,
                                      mesh.ngsolvemesh(),
                                      definedon=mesh.boundaries
                                                             ("contact"))

        P = ngsolve.Integrate(ngsolve.grad(potential) *
                              ngsolve.Conj(self.__sigma *
                                           ngsolve.grad(potential)),
                              mesh.ngsolvemesh())
        return electric_field, V_contact, P, potential
