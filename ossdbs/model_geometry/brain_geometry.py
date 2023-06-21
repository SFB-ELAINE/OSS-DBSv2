from .bounding_box import BoundingBox
import netgen.occ
import numpy as np
import logging

_logger = logging.getLogger(__name__)


class BrainGeometry:
    """CAD model of the brain geometry

    Parameters
    ----------

    shape: str
        Choose between "Sphere", "Ellipsoid", "Box", and "Custom"
    bounding_box : BoundingBox
        Bounding box of the geometry, can be none for a custom geometry

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
    def __init__(self,
                 shape: str,
                 bounding_box: BoundingBox = None
                 ) -> None:
        self._bbox = bounding_box
        self._geometry = None
        self._shape = shape

    @property
    def bounding_box(self):
        return self._bbox

    def get_surface_names(self):
        surface_list = []
        for face in self.geometry.faces:
            if face.name not in surface_list:
                surface_list.append(face.name)
        return surface_list

    @property
    def geometry(self) -> netgen.occ.Solid:
        """Return OCC shape of the brain geometry

        Returns
        -------
        netgen.occ.Solid
        """

        if self._geometry is None:
            _logger.debug("Inititalize geo of shape {}")
            self._geometry = self._create_shape()
            self._geometry.bc('BrainSurface')
            self._geometry.mat('Brain')

        return self._geometry

    def _create_shape(self) -> netgen.occ.Solid:
        if self._shape == "Sphere":
            return self._create_sphere()
        elif self._shape == "Box":
            return self._create_box()
        elif self._shape == "Ellipsoid":
            return self._create_ellipsoid()
        elif self._shape == "Custom":
            raise RuntimeError("Please load the custom geometry from a CAD file using the `import_geometry` or add it via the `set_geometry` method.")
        else:
            raise NotImplementedError("The shape `{}` is not implemented. Please choose among these shapes: Box, Sphere, Ellipsoid and Custom.".format(self._shape))

    def _create_ellipsoid(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self._bbox.end, self._bbox.start) / 2
        transformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return transformator(sphere).Move(self._bbox.start)

    def _create_sphere(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self._bbox.end, self._bbox.start) / 2
        center = np.add(self._bbox.start, (x, y, z))
        radius = np.min([x, y, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(center), r=radius)
        return sphere

    def _create_box(self) -> netgen.occ.Solid:
        box = netgen.occ.Box(self._bbox.start, self._bbox.end)
        return box

    def import_geometry(self, path_to_geo_file: str):
        """Import brain geometry from CAD file

        Notes
        -----

        TODO link to NGSolve / Netgen docs.

        """
        _logger.debug("Import brain geometry from file: {}".format(path_to_geo_file))
        occgeo = netgen.occ.OCCGeometry(path_to_geo_file)
        self._geometry = occgeo.shape
        bbox = self._geometry.bounding_box
        print(bbox)
        self._bbox = BoundingBox(bbox[0], bbox[1])
        self._geometry.bc('BrainSurface')
        self._geometry.mat('Brain')

    def set_geometry(self, geo: netgen.occ.Solid):
        """Set brain geometry from externally prepared OCC solid

        Notes
        -----


        This function is the only possibility to use custom-defined
        names for surface parts of the brain!
        """
        self._geometry = geo
        bbox = self._geometry.bounding_box
        self._bbox = BoundingBox(bbox[0], bbox[1])
        self._geometry.mat('Brain')
        # TODO test naming
        for face in self._geometry.faces:
            if face.name == "":
                face.name = "BrainSurface"
