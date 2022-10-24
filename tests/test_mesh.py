from src.mesh import Mesh
from src.volume_conductor_model import VolumeConductor
import netgen.occ as occ
import numpy as np

class GeometryDummy:
    def __init__(self) -> None:
        model = occ.Box(p1=occ.Pnt(0,0,0), p2=occ.Pnt(1,1,1))
        model.bc('contact')
        model.mat('saline')
        self.__geometry = occ.OCCGeometry(model)

    def ng_mesh(self):
        return self.__geometry.GenerateMesh()

def test_element_sizes():
    geo = GeometryDummy()
    mesh = Mesh(geometry=geo, order=2, boundaries={})
    np.testing.assert_allclose(mesh.element_sizes(), [1/6] * 12, atol=0.001)

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

def test_elements_in_region():
    geo = GeometryDummy()
    mesh = Mesh(geometry=geo, order=2, boundaries={})
    expected = [False, False, True, True, False, False, 
                True, False, False, True, False, False]
    assert mesh.elements_in_region((0, 0, 1), 0.5) == expected

def test_foo():
    geo = GeometryDummy()
    mesh = Mesh(geometry=geo, order=2, boundaries={"contact": 1.0, "wire": 0.0})
    mesh.foo()
    model = VolumeConductor(mesh, {"saline" : 1278*1e-6/1e-2,})
   
    model.evaluate_potential()

    field, contact, P, potential = model.evaluate_potential()
        
    
    print('field: ', field)
    print('voltage_contact: ', contact)
    print('impedance: ', 1 / P)