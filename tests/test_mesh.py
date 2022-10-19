from src.mesh import Mesh
import netgen.occ as occ
import numpy as np

class GeometryDummy:
    def __init__(self) -> None:
        model = occ.Box(p1=occ.Pnt(0,0,0), p2=occ.Pnt(1,1,1))
        self.__geometry = occ.OCCGeometry(model)

    def ng_mesh(self):
        return self.__geometry.GenerateMesh()

def test_foo():
    geo = GeometryDummy()
    mesh = Mesh(geo, 2, {})
    print(mesh.element_sizes())