
from src.electrodes.abstract_electrode import Electrode
from src.mesh import Mesh
from typing import List
import netgen
import numpy as np


class BrainModel:

    def __init__(self, region) -> None:
        self.__region = region
        self.__electrodes = []

    def add_electrodes(self, electrodes: List[Electrode]) -> None:
        for electrode in electrodes:
            self.__electrodes.append(electrode)

    def generate_mesh(self, order: int = 2) -> Mesh:
        geometry = self.__create_ellipsoid()
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()
        return Mesh(netgen.occ.OCCGeometry(geometry), order=order)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = trasformator(sphere).Move(self.__region.start)
        ellipsoid.bc('Brain')
        return ellipsoid
