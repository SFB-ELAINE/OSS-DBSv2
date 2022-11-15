from ossdbs.dielectric_model import DielectricModel1
from ossdbs.brainsubstance import Material
from ossdbs.electrodes import AbstractElectrode
from ossdbs.brain_imaging import MagneticResonanceImage
from ossdbs.brain_imaging import DiffusionTensorImage
from ossdbs.geometry import Geometry
from ossdbs.mesh import Mesh
import netgen
import numpy as np

from .voxel_space import VoxelSpace


class BrainModel:

    def __init__(self,
                 mri: MagneticResonanceImage,
                 dti: DiffusionTensorImage = None,
                 electrodes: list = None) -> None:
        self.__mri = mri
        self.__dti = dti
        self.__electrodes = electrodes if electrodes else []

    def bounding_box(self) -> tuple:
        return self.__mri.bounding_box()

    def material_distribution(self, material: Material) -> VoxelSpace:
        start, end = self.__mri.bounding_box()
        data = self.__material_distribution(material)
        return VoxelSpace(data=data, start=tuple(start), end=tuple(end))

    def __material_distribution(self, material: Material) -> np.ndarray:
        return self.__mri.data_map() == material

    def conductivity(self, frequency: float) -> None:
        csf_position = self.__material_distribution(Material.CSF)
        gm_position = self.__material_distribution(Material.GRAY_MATTER)
        wm_positiom = self.__material_distribution(Material.WHITE_MATTER)
        unknown_position = self.__material_distribution(Material.UNKNOWN)

        csf_model = DielectricModel1.create_model(Material.CSF)
        gm_model = DielectricModel1.create_model(Material.GRAY_MATTER)
        wm_model = DielectricModel1.create_model(Material.WHITE_MATTER)

        data = np.zeros(self.__mri.data_map().shape)
        data[csf_position] = csf_model.conductivity(frequency)
        data[gm_position] = gm_model.conductivity(frequency)
        data[wm_positiom] = wm_model.conductivity(frequency)
        data[unknown_position] = gm_model.conductivity(frequency)
        start, end = self.__mri.bounding_box()
        return VoxelSpace(data=data, start=tuple(start), end=tuple(end))

    def generate_mesh(self, order: int = 2):
        return Mesh(self.generate_geometry(), order=order)

    def generate_geometry(self):
        start, end = self.__mri.bounding_box()
        brain = Ellipsoid(start=start, end=end).create()

        geometry = brain
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()

        return BrainGeometry(geometry)


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


class BrainGeometry(Geometry):

    def __init__(self, geometry) -> None:
        self.__geometry = geometry

    def generate_mesh(self):
        return netgen.occ.OCCGeometry(self.__geometry).GenerateMesh()
