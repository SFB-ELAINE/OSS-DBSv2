
from ossdbs.electrodes.electrode import Electrode
from ossdbs.region import Region
from ossdbs.mesh import Mesh
from typing import List
import netgen
import numpy as np


class BrainGeometry:
    """Represents a simplification of the brain form.

    region : Region
        Location in 3D space of the geometry
    """
    def __init__(self, region: Region) -> None:
        self.__region = region
        self.__electrodes = []

    def set_electrodes(self, electrodes: List[Electrode]) -> None:
        """Set the electrodes, which are implemented in the brain.

        electrodes : list of Electrode objects
            Collection of electrodes.
        """
        self.__electrodes = [electrode for electrode in electrodes]

    def generate_mesh(self, order: int = 2) -> Mesh:
        """Generate a mesh based on the geometry and given mesh element order.

        order : int
            Order of mesh elements.
        """
        geometry = self.__create_ellipsoid()
        geometry.bc('Brain')
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()
        return Mesh(netgen.occ.OCCGeometry(geometry), order=order)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__region.start)
