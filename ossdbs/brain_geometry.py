
from ossdbs.electrodes.electrode import Electrode
from ossdbs.region import Region
from ossdbs.mesh import Mesh
from typing import List
import netgen
import ngsolve
import numpy as np


class BrainGeometry:
    """Represents a simplification of the brain form.

    region : Region
        Location in 3D space of the geometry
    """
    def __init__(self, region: Region, electrodes: List[Electrode]) -> None:
        self.__region = region
        self.__electrodes = electrodes

    def generate_mesh(self, parameters: dict = None, order: int = 2) -> Mesh:
        """Generate a mesh based on the geometry and given mesh element order.

        order : int
            Order of mesh elements.
        """
        netgen_geometry = self.__generate_geometry()
        ngmesh = netgen_geometry.GenerateMesh(parameters)
        mesh = ngsolve.Mesh(ngmesh=ngmesh)
        return Mesh(mesh, order=order)

    def load_mesh(self, file_name: str, order: int = 2) -> Mesh:
        netgen_geometry = self.__generate_geometry()
        mesh = ngsolve.Mesh(filename=file_name)
        mesh.SetGeometry(netgen_geometry)
        return Mesh(mesh, order=order)

    def __generate_geometry(self):
        geometry = self.__create_ellipsoid()
        geometry.bc('Brain')
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()

        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        return netgen_geometry

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__region.start)
