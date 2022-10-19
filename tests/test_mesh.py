from src.mesh import Mesh
import netgen.occ as occ
import numpy as np

class GeometryDummy:
    def __init__(self) -> None:
        model = occ.Box(p1=occ.Pnt(0,0,0), p2=occ.Pnt(1,1,1))
        self.__geometry = occ.OCCGeometry(model)

    def ng_mesh(self):
        return self.__geometry.GenerateMesh()

def test_element_sizes():
    geo = GeometryDummy()
    mesh = Mesh(geometry=geo, order=2, boundaries={})
    np.testing.assert_allclose(mesh.element_sizes(), [1/6]* 12, atol=0.001)

def test_centroids_of_elements():
    geo = GeometryDummy()
    mesh = Mesh(geometry=geo, order=2, boundaries={})
    expected = [(0.125, 0.625, 0.375),
                (0.875, 0.375, 0.375),
                (0.125, 0.375, 0.625),
                (0.375, 0.125, 0.375),
                (0.625, 0.875, 0.375),
                (0.375, 0.875, 0.625),
                (0.375, 0.375, 0.875),
                (0.625, 0.625, 0.875),
                (0.875, 0.625, 0.625),
                (0.625, 0.125, 0.625),
                (0.625, 0.375, 0.125),
                (0.375, 0.625, 0.125)] 

    np.testing.assert_allclose(mesh.centroids_of_elements(), 
                                expected,
                                atol=0.001)
