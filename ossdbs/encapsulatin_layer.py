
from ossdbs.bounding_box import BoundingBox
from ossdbs.electrode_models.electrode import ElectrodeModel
from ossdbs.mesh import Mesh
from typing import List
import numpy as np
import netgen
import ngsolve


class EncapsulatingLayers:

    def __init__(self,
                 electrodes: List[ElectrodeModel],
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

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        return sum([electrode.capsule_geometry(self.__thickness, self.__max_h)
                    for electrode in self.__electrodes])

    def bounding_boxes(self) -> List[BoundingBox]:
        return [self.__bounding_box(elec) for elec in self.__electrodes]

    def __bounding_box(self, electrode: ElectrodeModel) -> BoundingBox:
        geometry = electrode.capsule_geometry(self.__thickness)
        start, end = geometry.bounding_box
        return BoundingBox(tuple(start), tuple(end))

    def is_included(self, points: np.ndarray) -> np.ndarray:
        netgen_geometry = netgen.occ.OCCGeometry(self.geometry())
        ng_mesh = netgen_geometry.GenerateMesh()
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=2)
        return mesh.is_included(points)
