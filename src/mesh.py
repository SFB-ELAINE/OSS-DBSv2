from src.geometry import SimpleGeometry
import ngsolve


class Mesh:

    def __init__(self,
                geometry: SimpleGeometry, 
                order: int, 
                boundaries: dict) -> None:
        self.__mesh = ngsolve.Mesh(geometry.generate_mesh())
        self.__mesh.Curve(order)
        self.__order = order
        self.__boundaries = boundaries

    def ngsolvemesh(self):
        return self.__mesh

    def boundaries(self, name):
        return self.__mesh.Boundaries(name)

    def boundary_values(self) -> ngsolve.fem.CoefficientFunction:
        return self.__mesh.BoundaryCF(self.__boundaries)
        
    def get_materials(self) -> tuple:
        return self.__mesh.GetMaterials()

    def refine(self, error) -> None:
        errors = ngsolve.Integrate(error, self.__mesh, ngsolve.VOL,
                                    element_wise=True).real
        limit = 0.5 * max(errors)
        for element in self.__mesh.Elements():
            self.__mesh.SetRefinementFlag(element, errors[element.nr] > limit)
        self.__mesh.Refine()
        self.__mesh.Curve(self.__order)

    def space(self):
        dirichlet = '|'.join(str(key) for key in self.__boundaries.keys())
        return ngsolve.H1(self.__mesh, 
                            order=self.__order, 
                            dirichlet=dirichlet, 
                            complex=False, 
                            wb_withedges=False)

    def flux_space(self):
        return ngsolve.HDiv(self.__mesh, order=self.__order-1, complex=True) 

    def curve(self, order):
        self.__mesh.Curve(order)