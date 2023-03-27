
from ossdbs.electrodes.electrodes import Electrodes
from ossdbs.electrodes.encapsulation_layer import EncapsulatingLayers
from ossdbs.bounding_box import BoundingBox
import netgen
import ngsolve
import numpy as np

from ossdbs.mesh import Mesh


class BrainGeometry:
    """Represents a simplification of the brain form.

    bounding_box : BoundingBox
        Localization in 3D space of the geometry

    electrodes : Electrodes
        Collection of electode models

    encapsulation : EncapsulatingLayers
        Encapsulation of electrodes.
    """
    def __init__(self,
                 bounding_box: BoundingBox,
                 electrodes: Electrodes,
                 encapsulation: EncapsulatingLayers) -> None:
        self.__bbox = bounding_box
        self.__electrodes = electrodes
        self.__encapsulating_layer = encapsulation

    def geometry(self) -> netgen.libngpy._NgOCC.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.libngpy._NgOCC.OCCGeometry
        """
        body = self.__create_ellipsoid()
        body.bc('BrainSurface')
        electrodes_geometry = self.__electrodes.geometry()

        if not self.__encapsulating_layer.thickness:
            return netgen.occ.OCCGeometry(body - electrodes_geometry)

        capsule_geometry = self.__encapsulating_layer.geometry()
        cut = capsule_geometry + electrodes_geometry - body
        part_1 = body - capsule_geometry - electrodes_geometry
        part_2 = capsule_geometry - electrodes_geometry - cut
        geometry = netgen.occ.Glue([part_1, part_2])
        return netgen.occ.OCCGeometry(geometry)

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
        netgen_geometry = netgen.occ.OCCGeometry(self.geometry())
        ng_mesh = netgen_geometry.GenerateMesh()
        ngsolve_mesh = ngsolve.Mesh(ngmesh=ng_mesh)
        mesh = Mesh(ngsolve_mesh=ngsolve_mesh, order=2)
        return mesh.is_included(points)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__bbox.end, self.__bbox.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__bbox.start)
