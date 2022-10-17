from src.geometry import SimpleGeometry
import ngsolve


class Mesh:

    def __init__(self,
                geometry: SimpleGeometry, 
                order: int, 
                boundaries: dict) -> None:
        self.__mesh = ngsolve.Mesh(ngmesh=geometry.generate_mesh())
        self.__mesh.Curve(order=order)
        self.__order = order
        self.__boundaries = boundaries

    def boundaries(self, name: str) -> ngsolve.comp.Region:
        return self.__mesh.Boundaries(pattern=name)

    def boundary_coefficients(self) -> ngsolve.fem.CoefficientFunction:
        return self.__mesh.BoundaryCF(values=self.__boundaries)

    def curve(self, order: int) -> None:
        self.__mesh.Curve(order=order)

    def flux_space(self) -> ngsolve.comp.HDiv:
        return ngsolve.HDiv(mesh=self.__mesh, 
                            order=self.__order-1,
                            complex=True)    

    def materials(self) -> tuple:
        return self.__mesh.GetMaterials()

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        return self.__mesh

    def refine(self, error: ngsolve.fem.CoefficientFunction) -> None:
        errors = ngsolve.Integrate(cf=error,
                                    mesh=self.__mesh, 
                                    VOL_or_BND=ngsolve.VOL,
                                    element_wise=True
                                    ).real
        limit = 0.5 * max(errors)
        for element in self.__mesh.Elements():
            self.__mesh.SetRefinementFlag(ei=element, 
                                            refine=errors[element.nr] > limit)
        
        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def space(self) -> ngsolve.comp.H1:
        dirichlet = '|'.join(str(key) for key in self.__boundaries.keys())
        return ngsolve.H1(mesh=self.__mesh, 
                            order=self.__order, 
                            dirichlet=dirichlet, 
                            complex=False, 
                            wb_withedges=False)
                            