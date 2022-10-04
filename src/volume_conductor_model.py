
from src.mesh import Mesh
from src.laplace_equation import LaplaceEquation
import ngsolve


class VolumeConductor:
    
    def __init__(self, mesh: Mesh, conductivity: dict) -> None:
        self.__mesh = mesh
        self.__sigma = ngsolve.CoefficientFunction([conductivity[mat] 
                                             for mat in mesh.materials()]) 
        
    def evaluate_potential(self, 
                            ndof_limit : int = 1e6, 
                            n_refinements: int  = 1) -> tuple:
        potential = self.__solve_bvp()
        iterations = 0
        while (self.__mesh.space().ndof < ndof_limit and 
                                                iterations < n_refinements):
            self.__mesh.refine(self.__error(potential))
            potential = self.__solve_bvp()
            iterations = iterations + 1
        return self.__postprocess(potential)

    def __solve_bvp(self) -> ngsolve.comp.GridFunction:
        boundaries = self.__mesh.boundary_coefficients()
        potential = ngsolve.GridFunction(self.__mesh.space())
        potential.Set(boundaries, ngsolve.BND)
        equation = LaplaceEquation(self.__mesh.space(),  self.__sigma)
        potential.vec.data = equation.solve_bvp(potential)
        return potential

    def __error(self, potential: ngsolve.comp.GridFunction) \
                                            -> ngsolve.fem.CoefficientFunction:
        flux = ngsolve.grad(potential)
        flux_potential = ngsolve.GridFunction(self.__mesh.flux_space())
        flux_potential.Set(flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)

    def __postprocess(self, potential: ngsolve.comp.GridFunction) -> tuple:
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

   