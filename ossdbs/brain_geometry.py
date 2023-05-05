from ossdbs.bounding_box import BoundingBox
import netgen
import numpy as np
import logging

_logger = logging.getLogger(__name__)


class BrainGeometry:
    """CAD model of the brain geometry

    Parameters
    ----------

    bounding_box : BoundingBox
        Bounding box of the geometry
    shape: str
        Choose between "Sphere", "Ellipsoid", "Box", and "Custom"

    Notes
    -----

    The geometry will be initialized once.
    A custom shape is initialized with the argument `Custom` but
    requires to load the geometry from a STEP or BREP file using
    `import_geometry`.


    The brain has a single surface named `BrainSurface` and
    a single volume named `Brain`.


    Examples
    --------

    TODO example for custom import

    """
    def __init__(self,
                 bounding_box: BoundingBox,
                 shape: str
                 ) -> None:
        self.__bbox = bounding_box
        self.__geometry = None
        self.__shape = shape

    @property
    def bounding_box(self):
        return self.__bbox

    @property
    def geometry(self) -> netgen.occ.Solid:
        """Return OCC shape of the brain geometry

        Returns
        -------
        netgen.occ.Solid
        """

        if self.__geometry is None:
            _logger.debug("Inititalize geo of shape {}")
            self.__geometry = self.__create_shape()
            self.__geometry.bc('BrainSurface')
            self.__geometry.mat('Brain')

        return self.__geometry

    def __create_shape(self) -> netgen.occ.Solid:
        if self.__shape == "Sphere":
            return self.__create_sphere()
        elif self.__shape == "Box":
            return self.__create_box()
        elif self.__shape == "Ellipsoid":
            return self.__create_ellipsoid()
        elif self.__shape == "Custom":
            raise RuntimeError("Please load the custom geometry from a CAD file using the `import_geometry` method.")
        else:
            raise NotImplementedError("The shape `{}` is not implemented. Please choose among these shapes: Box, Sphere, Ellipsoid and Custom.".format(self.__shape))

    def __create_ellipsoid(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self.__bbox.end, self.__bbox.start) / 2
        transformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return transformator(sphere).Move(self.__bbox.start)

    def __create_sphere(self) -> netgen.occ.Solid:
        x, y, z = np.subtract(self.__bbox.end, self.__bbox.start) / 2
        center = np.add(self.__bbox.start, (x, y, z))
        radius = np.min([x, y, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(center), r=radius)
        return sphere

    def __create_box(self) -> netgen.occ.Solid:
        box = netgen.occ.Box(self.__bbox.start, self.__bbox.end)
        return box

    def import_geometry(self, path_to_geo_file: str):
        """Import brain geometry from CAD file

        Notes
        -----

        TODO link to NGSolve / Netgen docs.

        """
        _logger.debug("Import brain geometry from file: {}".format(path_to_geo_file))
        occgeo = netgen.occ.OCCGeometry(path_to_geo_file)
        self.__geometry = occgeo.shape
        bbox = self.__geometry.bounding_box
        print(bbox)
        self.__bbox = BoundingBox(bbox[0], bbox[1])
        self.__geometry.bc('BrainSurface')
        self.__geometry.mat('Brain')

    def set_geometry(self, geo: netgen.occ.Solid):
        """Set brain geometry from externally prepared OCC solid
        """
        self.__geometry = geo
        bbox = self.__geometry.bounding_box
        self.__bbox = BoundingBox(bbox[0], bbox[1])
        self.__geometry.bc('BrainSurface')
        self.__geometry.mat('Brain')
