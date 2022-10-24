from xml.dom.minidom import Element
from src.geometry import SimpleGeometry
import ngsolve
import numpy as np


class Mesh:

    def __init__(self,
                geometry: SimpleGeometry, 
                order: int, 
                boundaries: dict) -> None:
        self.__mesh = ngsolve.Mesh(ngmesh=geometry.ng_mesh())
        self.__mesh.Curve(order=order)
        self.__order = order
        self.__boundaries = boundaries

    def boundaries(self, name: str) -> ngsolve.comp.Region:
        return self.__mesh.Boundaries(pattern=name)

    def boundary_coefficients(self) -> ngsolve.fem.CoefficientFunction:
        return self.__mesh.BoundaryCF(values=self.__boundaries)

    def centroids_of_elements(self) -> list:
        shape = (self.__mesh.ne, 4, 3)
        vertices = np.array([self.__mesh[v].point 
                             for element in self.__mesh.Elements()
                             for v in element.vertices]).reshape(shape)
        return [list(c) for c in np.sum(vertices, axis=1) / 4]

    def flux_space(self) -> ngsolve.comp.HDiv:
        return ngsolve.HDiv(mesh=self.__mesh, 
                            order=self.__order-1,
                            complex=True)    

    def materials(self) -> tuple:
        return self.__mesh.GetMaterials()

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        return self.__mesh

    def refine(self) -> None:
        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def refine_by_error(self, error: ngsolve.fem.CoefficientFunction) -> None:
        self.mark_elements_by_error(error)
        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def mark_elements_by_error(self, 
                               error: ngsolve.fem.CoefficientFunction) -> None:
        errors = ngsolve.Integrate(cf=error,
                                    mesh=self.__mesh, 
                                    VOL_or_BND=ngsolve.VOL,
                                    element_wise=True).real
        limit = 0.5 * max(errors)
        flags = [errors[el.nr] > limit for el in self.__mesh.Elements()] 
        for index, element in enumerate(self.__mesh.Elements()):
            self.__mesh.SetRefinementFlag(ei=element, refine=flags[index])

    def element_sizes(self) -> list:
        volumes = ngsolve.Integrate(cf=ngsolve.CoefficientFunction(1), 
                                    mesh=self.__mesh,
                                    element_wise=True).NumPy()
        return list((6 * volumes) ** 1 / 3)
        
    def space(self) -> ngsolve.comp.H1:
        dirichlet = '|'.join(str(key) for key in self.__boundaries.keys())
        return ngsolve.H1(mesh=self.__mesh, 
                          order=self.__order, 
                          dirichlet=dirichlet, 
                          complex=False, 
                          wb_withedges=False)

    def elements_in_region(self, point: tuple, r: int) -> np.ndarray:
        shape = (self.__mesh.ne, 4, 3)
        vertices = np.array([self.__mesh[v].point 
                             for element in self.__mesh.Elements()
                             for v in element.vertices]).reshape(shape)

        distance = np.sum((vertices - point) ** 2, axis=2)
        return list(np.any(distance <= r ** 2, axis=1))

    def foo(self):

        # help( ngsolve.VoxelCoefficient)
        data = np.array([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]).reshape(2,2,2)
        #data = np.array([np.eye(3)]*8).reshape(2,2,2,3,3)
        
        cf = ngsolve.VoxelCoefficient(start=(0,0,0), end=(1,1,1), values=data, linear=False)

        volumes = ngsolve.Integrate(cf=cf, 
                                    mesh=self.__mesh,
                                    element_wise=True).NumPy()
        print(volumes)
 

