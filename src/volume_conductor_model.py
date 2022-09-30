
from src.mesh import Mesh
import ngsolve


class VolumeConductor:
    
    def __init__(self, mesh: Mesh, conductivity) -> None:
        self.__mesh = mesh
        self.__sigma = ngsolve.CoefficientFunction([conductivity[mat] 
                                             for mat in mesh.get_materials()]) 
        
    def evaluate_potential(self, ndof_limit = 1e6, n_refinements = 1):
        potential = self.__solve_bvp()
        iterations = 0
        while (self.__mesh.space().ndof < ndof_limit and 
                                                iterations < n_refinements):
            self.__mesh.refine(self.__error(potential))
            potential = self.__solve_bvp()
            iterations = iterations + 1
        return self.__postprocess(potential)

    def __solve_bvp(self):
        boundaries = self.__mesh.boundary_values()
        potential = ngsolve.GridFunction(self.__mesh.space())
        potential.Set(boundaries, ngsolve.BND)
        equation = LaplaceEquation(self.__mesh.space(),  self.__sigma)
        potential.vec.data = equation.solve_bvp(potential)
        return potential

    def __error(self, potential) -> list:
        flux = ngsolve.grad(potential)
        flux_potential = ngsolve.GridFunction(self.__mesh.flux_space())
        flux_potential.Set(flux)
        difference = flux - flux_potential
        return difference * ngsolve.Conj(difference)

    def __postprocess(self, potential):
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

    
class LaplaceEquation:

    def __init__(self, space, factor) -> None:
        u = space.TrialFunction()
        v = space.TestFunction()
        self.__a = ngsolve.BilinearForm(space, symmetric=True)
        self.__a += factor * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        self.__f = ngsolve.LinearForm(space)
        self.__preconditioner = ngsolve.Preconditioner(self.__a, 
                                                        type="bddc",
                                                        coarsetype="h1amg"
                                                        )

    def solve_bvp(self, input):
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