
from ossdbs.bounding_box import BoundingBox
from ossdbs.mesh import Mesh
import numpy as np
import netgen
import ngsolve


class EncapsulatingLayers:

    def __init__(self,
                 electrodes: list,
                 thickness: float
                 ) -> None:
        self.__electrodes = electrodes
        self.__thickness = thickness
        self.__max_h = 0.1

    @property
    def thickness(self) -> float:
        return self.__thickness

    def set_max_h(self, max_h: float = 0.1) -> None:
        self.__max_h = max_h

    def geometry(self):
        return sum([electrode.capsule_geometry(self.__thickness, self.__max_h)
                    for electrode in self.__electrodes])

    def bounding_boxes(self):
        bounding_boxes = []
        for electrode in self.__electrodes:
            geometry = electrode.capsule_geometry(self.__thickness)
            start, end = geometry.bounding_box
            bounding_boxes.append(BoundingBox(tuple(start), tuple(end)))
        return bounding_boxes

    def is_included(self, points: np.ndarray) -> np.ndarray:
        geometry = self.geometry()
        netgen_geometry = netgen.occ.OCCGeometry(geometry)
        ng_mesh = netgen_geometry.GenerateMesh()
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=2)
        return mesh.is_included(points)
