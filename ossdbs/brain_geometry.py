
from ossdbs.electrodes.electrode import Electrode
from ossdbs.region import BoundingBox
from typing import List
import netgen
import numpy as np


class BrainGeometry:
    """Represents a simplification of the brain form.

    region : Region
        Location in 3D space of the geometry
    """
    def __init__(self, region: BoundingBox, electrodes: List[Electrode])\
            -> None:
        self.__region = region
        self.__electrodes = electrodes

    def netgen_geometry(self) -> netgen.libngpy._NgOCC.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.libngpy._NgOCC.OCCGeometry
        """

        geometry = self.__create_ellipsoid()
        #  geometry = self.__create_box()
        geometry.bc('BrainSurface')
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()

        for electrode in self.__electrodes:
            capsule_geo = electrode.encapsulating_geometry(thickness=0.0001)
            electrode_geo = electrode.generate_geometry()
            geometry = netgen.occ.Glue([geometry - capsule_geo,
                                        capsule_geo - electrode_geo])

        return netgen.occ.OCCGeometry(geometry)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__region.start)

    def __create_box(self) -> netgen.libngpy._NgOCC.Solid:
        return netgen.occ.Box(self.__region.start, self.__region.end)
