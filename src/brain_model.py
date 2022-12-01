from src.brain_imaging import MagneticResonanceImage
from src.brain_imaging import DiffusionTensorImage
from src.geometry import Geometry
from src.mesh import Mesh
import netgen
import numpy as np


class BrainModel:

    def __init__(self,
                 mri: MagneticResonanceImage,
                 dti: DiffusionTensorImage = None,
                 electrodes: list = None) -> None:
        self.__mri = mri
        self.__electrodes = electrodes if electrodes else []

    def generate_mesh(self, order: int = 2):
        return Mesh(self.generate_geometry(), order=order)

    def generate_geometry(self):
        start, end = self.__mri.bounding_box()
        brain = Ellipsoid(start=start, end=end).create()

        geometry = brain
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()

        return Geometry(geometry=geometry)


class Ellipsoid:

    def __init__(self, start: tuple, end: tuple) -> None:
        self.__start = start
        self.__end = end

    def create(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        x, y, z = (np.array(self.__end) - np.array(self.__start)) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = trasformator(sphere).Move(tuple(self.__start))
        ellipsoid.bc('Brain')
        return ellipsoid
