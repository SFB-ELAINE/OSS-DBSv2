
import netgen


class Geometry():

    def __init__(self, geometry) -> None:
        self.__geometry = geometry

    def generate_mesh(self):
        return netgen.occ.OCCGeometry(self.__geometry).GenerateMesh()
