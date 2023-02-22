
from ossdbs.electrodes.electrode import Electrode
from ossdbs.region import BoundingBox
from typing import List
import netgen
import netgen.occ as occ
import numpy as np


class BrainGeometry:
    """Represents a simplification of the brain form.

    region : Region
        Location in 3D space of the geometry
    """
    def __init__(self,
                 region: BoundingBox,
                 electrodes: List[Electrode],
                 encapsulating_thickness: float,
                 ) -> None:
        self.__region = region
        self.__electrodes = electrodes
        self.__encapsulating_thickness = encapsulating_thickness
        self.__edges = []
        self.__max_h = 0.5

    def set_edge_max_h(self, max_h: float) -> None:
        self.__max_h = max_h

    def edges_for_finer_meshing(self, edges: List[str]) -> None:
        self.__edges = edges

    def netgen_geometry(self) -> netgen.libngpy._NgOCC.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.libngpy._NgOCC.OCCGeometry
        """
        body = self.__create_ellipsoid()
        #  geometry = self.__create_box()
        body.bc('BrainSurface')
        thickness = self.__encapsulating_thickness

        elec_geo = sum([elec.geometry() for elec in self.__electrodes])

        if not thickness:
            geometry = body - elec_geo
            self.__set_maxh_edges(geometry)
            return occ.OCCGeometry(geometry)

        capsule_geo = sum([electrode.capsule_geometry(thickness)]
                          for electrode in self.__electrodes)
        cut = capsule_geo - body
        geometry = occ.Glue([body - capsule_geo, capsule_geo - elec_geo - cut])
        self.__set_maxh_edges(geometry)
        return netgen.occ.OCCGeometry(geometry)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__region.start)

    def __create_box(self) -> netgen.libngpy._NgOCC.Solid:
        return netgen.occ.Box(self.__region.start, self.__region.end)

    def __set_maxh_edges(self,
                         geometry: netgen.libngpy._NgOCC.Solid,
                         ) -> None:
        for edge in geometry.edges:
            if edge.name in self.__edges:
                edge.maxh = self.__max_h
