
from ossdbs.bounding_box import BoundingBox
from .electrode_models import ElectrodeModel
from ossdbs.fem.mesh import Mesh
from typing import List
import numpy as np
import netgen
import ngsolve


class EncapsulatingLayers:
    """Collection of electrode encapsulating layers.

    Attributes
    ----------
    electrodes : list of ElectrodeModel
        Collection of electrode models.

    thickness : float
        Thickness of encapsulating layers around electrodes.
    """

    def __init__(self,
                 electrode_models: List[ElectrodeModel],
                 thickness: float
                 ) -> None:
        self.__electrodes = electrode_models
        self.__thickness = thickness
        self.__max_h = 0.1

    @property
    def thickness(self) -> float:
        """Return the thickness of the encapsulating layers."""
        return self.__thickness

    def bounding_boxes(self) -> List[BoundingBox]:
        """Return a collection of bounding boxes. A bounding box for each
        enapsulation layer.

        Returns
        -------
        list of Bounding Boxes
        """
        return [self.__bounding_box(elec) for elec in self.__electrodes]

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """"Return the geometry of the combined encapsulating layers.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        geometry = sum([electrode.capsule_geometry(self.__thickness)
                        for electrode in self.__electrodes])
        geometry.maxh = self.__max_h
        return geometry

    def is_included(self, points: np.ndarray) -> np.ndarray:
        """Check each point in collection for collision with geometry.
        True if point is included in geometry, false otherwise.

        Parameters
        ----------
        points: np.ndarray
            Array of point coordinates (x, y, z).

        Returns
        -------
        np.ndarray
            Array representing the state of collision for each point.
            True if point is included in geometry, False otherwise.
        """
        if not self.thickness:
            return np.full(len(points), False)

        netgen_geometry = netgen.occ.OCCGeometry(self.geometry())
        ng_mesh = netgen_geometry.GenerateMesh()
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=2)
        return mesh.is_included(points)

    def set_max_h(self, max_h: float = 0.1) -> None:
        """Set the maximum height for a mesh element in these encapsulating
        layers.

        Parameters
        ----------
        max_h : float
            Maximum height for a mesh element.
        """
        self.__max_h = max_h

    def __bounding_box(self, electrode: ElectrodeModel) -> BoundingBox:
        geometry = electrode.capsule_geometry(self.__thickness)
        start, end = geometry.bounding_box
        return BoundingBox(tuple(start), tuple(end))
