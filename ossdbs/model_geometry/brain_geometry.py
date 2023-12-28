import logging

import netgen.occ
import numpy as np
import numpy.linalg as npl

from .bounding_box import BoundingBox

_logger = logging.getLogger(__name__)


class BrainGeometry:
    """CAD model of the brain geometry.

    Parameters
    ----------
    shape: str
        Choose between "Sphere", "Ellipsoid", "Box", and "Custom"
    bounding_box : BoundingBox
        Bounding box of the geometry in voxel space, can be none for custom geometry

    Notes
    -----
    The geometry will be initialized once.
    A custom shape can be loaded from a STEP or BREP file using
    `import_geometry`.

    The brain has a single surface named `BrainSurface` and
    a single volume named `Brain`.
    An exception can be made by adding the geometry via the
    `set_geometry` method.


    Examples
    --------
    TODO example for custom import

    """

    def __init__(
        self,
        shape: str,
        bounding_box: BoundingBox = None,
        trafo_matrix: np.ndarray = None,
        translation: np.ndarray = None,
    ) -> None:
        self._bbox = bounding_box
        self._geometry = None
        self._shape = shape
        if trafo_matrix is None:
            self._trafo_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        else:
            self._trafo_matrix = trafo_matrix
        if translation is None:
            self._translation = np.array([0, 0, 0])
        else:
            self._translation = translation
        self._affine_trafo = netgen.occ.gp_GTrsf(
            mat=self._trafo_matrix.ravel(), vec=self._translation
        )

    def get_surface_names(self):
        """Get names of surfaces in geometry."""
        surface_list = []
        for face in self.geometry.faces:
            if face.name not in surface_list:
                surface_list.append(face.name)
        return surface_list

    @property
    def geometry(self) -> netgen.occ.Solid:
        """Return OCC shape of the brain geometry.

        Returns
        -------
        netgen.occ.Solid
        """
        if self._geometry is None:
            _logger.debug(f"Inititalize geo of shape {self._shape}")
            self._geometry = self._create_shape()
            self._geometry.bc("BrainSurface")
            self._geometry.mat("Brain")

        return self._geometry

    def _create_shape(self) -> netgen.occ.Solid:
        if self._shape == "Sphere":
            return self._create_sphere()
        elif self._shape == "Box":
            return self._create_box()
        elif self._shape == "Ellipsoid":
            return self._create_ellipsoid()
        elif self._shape == "Custom":
            raise RuntimeError(
                """Please load the custom geometry from a CAD file using the
                `import_geometry` or add it via the `set_geometry` method."""
            )
        else:
            raise NotImplementedError(
                f"""The shape {self._shape} is not implemented.
                Please choose among these shapes: Box, Sphere, Ellipsoid and Custom."""
            )

    def _create_ellipsoid(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self._bbox.end, self._bbox.start) / 2
        transformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = transformator(sphere).Move(self._bbox.start)

        return self._affine_trafo(ellipsoid)

    def _create_sphere(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self._bbox.end, self._bbox.start) / 2
        center = np.add(self._bbox.start, (x, y, z))
        radius = np.min([x, y, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(center), r=radius)
        return self._affine_trafo(sphere)

    def _create_box(self) -> netgen.occ.Solid:
        box = netgen.occ.Box(self._bbox.start, self._bbox.end)
        return self._affine_trafo(box)

    def import_geometry(self, path_to_geo_file: str):
        """Import brain geometry from CAD file.

        Notes
        -----
        TODO link to NGSolve / Netgen docs.

        """
        _logger.debug(f"Import brain geometry from file: {path_to_geo_file}")
        occgeo = netgen.occ.OCCGeometry(path_to_geo_file)
        self._geometry = occgeo.shape
        self._geometry.bc("BrainSurface")
        self._geometry.mat("Brain")
        # get bounding box in voxel space
        bbox = self._geometry.bounding_box
        inv_trafo = npl.inv(self._trafo_matrix)
        inv_affine_trafo = netgen.occ.gp_GTrsf(
            mat=inv_trafo.ravel(), vec=-inv_trafo.dot(self._translation)
        )
        box = netgen.occ.Box(bbox[0], bbox[1])
        box_voxel = inv_affine_trafo(box)
        bbox_voxel = box_voxel.bounding_box
        self._bbox = BoundingBox(bbox_voxel[0], bbox_voxel[1])

    def set_geometry(self, geo: netgen.occ.Solid):
        """Set brain geometry from externally prepared OCC solid.

        Notes
        -----
        This function is the only possibility to use custom-defined
        names for surface parts of the brain!
        """
        self._geometry = geo
        bbox = self._geometry.bounding_box
        self._bbox = BoundingBox(bbox[0], bbox[1])
        self._geometry.mat("Brain")
        # TODO test naming
        for face in self._geometry.faces:
            if face.name == "":
                face.name = "BrainSurface"
