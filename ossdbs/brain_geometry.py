
from ossdbs.electrode_collection import Electrodes
from ossdbs.encapsulatin_layer import EncapsulatingLayers
from ossdbs.bounding_box import BoundingBox
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
                 electrodes: Electrodes,
                 encap_layer: EncapsulatingLayers) -> None:
        self.__region = region
        self.__electrodes = electrodes
        self.__encapsulating_layer = encap_layer

    def netgen_geometry(self) -> netgen.libngpy._NgOCC.OCCGeometry:
        """Create a netgen geometry of this brain model.

        Returns
        -------
        netgen.libngpy._NgOCC.OCCGeometry
        """
        body = self.__create_ellipsoid()
        body.bc('BrainSurface')
        elec_geo = self.__electrodes.geometry()

        if not self.__encapsulating_layer.thickness:
            return occ.OCCGeometry(body - elec_geo)

        capsule_geo = self.__encapsulating_layer.geometry()
        cut = capsule_geo - body
        geometry = occ.Glue([body - capsule_geo, capsule_geo - elec_geo - cut])
        return netgen.occ.OCCGeometry(geometry)

    def __create_ellipsoid(self) -> netgen.libngpy._NgOCC.Solid:
        x, y, z = np.subtract(self.__region.end, self.__region.start) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        return trasformator(sphere).Move(self.__region.start)
