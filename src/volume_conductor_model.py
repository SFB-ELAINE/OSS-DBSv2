
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
    
    def __init__(self, mesh: Mesh, conductivity: dict) -> None:
        self.__mesh = mesh
        conductivities = [conductivity[mat] for mat in mesh.materials()]
        self.__sigma = ngsolve.CoefficientFunction(coef=conductivities) 
        
    def evaluate_potential(self, 
                            n_dof_limit : int = 1e6, 
                            n_refinements: int  = 1) -> tuple:
        """Evaluate electrical potential of volume conductor.
        
        Parameters
        ----------
        n_dof_limit : int
            Maximal number of Demensions of Freedom of FEM space.

        n_refinements: int
            Number of mesh refinements.

        Returns
        -------

        return : tuple
            Postprocessed data: lectric_field, V_contact, Power, potential
        
        """
        potential = self.__solve_bvp()
        iterations = 0
        while (self.__mesh.space().ndof < n_dof_limit and 
                                                iterations < n_refinements):
            self.__mesh.refine(error=self.__error(potential))
            potential = self.__solve_bvp()
            iterations = iterations + 1
        return self.__postprocess(potential=potential)

    def __solve_bvp(self) -> ngsolve.comp.GridFunction:
        potential = ngsolve.GridFunction(space=self.__mesh.space())
        potential.Set(coefficient=self.__mesh.boundary_coefficients(),
                      VOL_or_BND=ngsolve.BND)
        equation = LaplaceEquation(space=self.__mesh.space(), 
                                   coefficient=self.__sigma)
        potential.vec.data = equation.solve_bvp(input=potential)
        return potential

    def __error(self, potential: ngsolve.comp.GridFunction) \
                                            -> ngsolve.fem.CoefficientFunction:
        flux = ngsolve.grad(potential)
        flux_potential = ngsolve.GridFunction(space=self.__mesh.flux_space())
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)

    def __postprocess(self, potential: ngsolve.comp.GridFunction) -> tuple:
        # MappedIntegrationPoint
        electric_field = ngsolve.sqrt(ngsolve.grad(potential) * \
                                                        ngsolve.grad(potential))
        electric_field = electric_field(
                            self.__mesh.ngsolvemesh()(0, 0, 0.02 / 2)) * 1e3  
        V_contact = ngsolve.Integrate(potential,
                                      self.__mesh.ngsolvemesh(), 
                                      definedon=self.__mesh.boundaries
                                                                    ("contact"))

        P = ngsolve.Integrate(ngsolve.grad(potential) * \
                        ngsolve.Conj(self.__sigma * ngsolve.grad(potential)), 
                        self.__mesh.ngsolvemesh())
        return electric_field, V_contact, P, potential
