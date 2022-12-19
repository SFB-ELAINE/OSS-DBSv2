from src.geometry import Geometry
import ngsolve
import numpy as np

from src.voxels import Voxels


class Mesh:

    def __init__(self,
                 geometry: Geometry,
                 order: int) -> None:
        self.__mesh = ngsolve.Mesh(ngmesh=geometry.generate_mesh())
        self.__mesh.Curve(order=order)
        self.__order = order

    def boundaries(self, name: str) -> ngsolve.comp.Region:
        return self.__mesh.Boundaries(pattern=name)

    def get_boundaries(self):
        return list(set(self.__mesh.GetBoundaries()) - set(['default']))

    def boundary_coefficients(self, boundaries) \
            -> ngsolve.fem.CoefficientFunction:
        return self.__mesh.BoundaryCF(values=boundaries)

    def flux_space(self, complex: bool = True) -> ngsolve.comp.HDiv:
        return ngsolve.HDiv(mesh=self.__mesh,
                            order=self.__order-1,
                            complex=complex)

    def materials(self) -> tuple:
        return self.__mesh.GetMaterials()

    def ngsolvemesh(self) -> ngsolve.comp.Mesh:
        return self.__mesh

    def refine(self) -> None:
        self.__mesh.Refine()
        self.__mesh.Curve(order=self.__order)

    def elements_with_error(self,
                            error: ngsolve.fem.CoefficientFunction) -> None:
        errors = ngsolve.Integrate(cf=error,
                                   mesh=self.__mesh,
                                   VOL_or_BND=ngsolve.VOL,
                                   element_wise=True).real
        limit = 0.5 * max(errors)
        return [errors[el.nr] > limit for el in self.__mesh.Elements()]

    def elements_at_position(self, position: Voxels) -> None:
        space = ngsolve.L2(self.__mesh, order=0)
        grid_function = ngsolve.GridFunction(space=space)
        cf = ngsolve.VoxelCoefficient(start=tuple(position.start),
                                      end=tuple(position.end),
                                      values=position.data.astype(float),
                                      linear=False)
        grid_function.Set(cf)
        return grid_function.vec.FV().NumPy()

    def elements_grater_than(self, limit: float) -> list:
        volumes = ngsolve.Integrate(cf=ngsolve.CoefficientFunction(1),
                                    mesh=self.__mesh,
                                    element_wise=True).NumPy()
        heights = (6 * volumes) ** 1 / 3
        return heights > limit

    def sobolev_space(self, complex: bool = False) -> ngsolve.comp.H1:
        dirichlet = '|'.join(boundary for boundary in self.get_boundaries())
        return ngsolve.H1(mesh=self.__mesh,
                          order=self.__order,
                          dirichlet=dirichlet,
                          complex=complex,
                          wb_withedges=False)

    def set_refinement_flag(self, flags):
        for index, element in enumerate(self.__mesh.Elements()):
            self.__mesh.SetRefinementFlag(ei=element, refine=flags[index])

    def set_boundary_refinement_flags(self, boundaries: list):
        for element in self.__mesh.Elements(ngsolve.BND):
            flag = element.mat in boundaries
            self.__mesh.SetRefinementFlag(ei=element, refine=flag)

    def centroids_of_elements(self) -> list:
        shape = (self.__mesh.ne, 4, 3)
        vertices = np.array([self.__mesh[v].point
                             for element in self.__mesh.Elements()
                             for v in element.vertices]).reshape(shape)
        return [list(c) for c in np.sum(vertices, axis=1) / 4]

    def element_sizes(self) -> list:
        volumes = ngsolve.Integrate(cf=ngsolve.CoefficientFunction(1),
                                    mesh=self.__mesh,
                                    element_wise=True).NumPy()
        return list((6 * volumes) ** 1 / 3)

    def bounding_box(self):
        points = [vertice.point for vertice in self.__mesh.vertices]
        start = tuple(np.min(points, axis=0))
        end = tuple(np.max(points, axis=0))
        return start, end
