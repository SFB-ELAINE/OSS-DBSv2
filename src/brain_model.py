from dataclasses import dataclass
from src.dielectric_model import DielectricModel1
from src.brainsubstance import Material
from src.brain_imaging import MagneticResonanceImage
from src.brain_imaging import DiffusionTensorImage
from src.geometry import Geometry
from src.mesh import Mesh
import netgen
import numpy as np
import ngsolve
from src.voxel_space import VoxelSpace


class BrainModel:

    def __init__(self,
                 mri: MagneticResonanceImage,
                 dti: DiffusionTensorImage = None,
                 electrodes: list = None) -> None:
        self.__mri = mri
        self.__electrodes = electrodes if electrodes else []

    def bounding_box(self) -> tuple:
        return self.__mri.bounding_box()

    def material_distribution(self, material: Material) -> VoxelSpace:
        start, end = self.__mri.bounding_box()
        data = self.__material_distribution(material)
        return VoxelSpace(data=data, start=tuple(start), end=tuple(end))

    def __material_distribution(self, material: Material) -> np.ndarray:
        return self.__mri.data_map() == material

    def conductivity(self, frequency: float) -> VoxelSpace:
        csf_position = self.__material_distribution(Material.CSF)
        gm_position = self.__material_distribution(Material.GRAY_MATTER)
        wm_position = self.__material_distribution(Material.WHITE_MATTER)

        csf_model = DielectricModel1.create_model(Material.CSF)
        gm_model = DielectricModel1.create_model(Material.GRAY_MATTER)
        wm_model = DielectricModel1.create_model(Material.WHITE_MATTER)

        default = gm_model.conductivity(frequency)
        data = np.full(self.__mri.data_map().shape, default)
        data[csf_position] = csf_model.conductivity(frequency)
        data[gm_position] = gm_model.conductivity(frequency)
        data[wm_position] = wm_model.conductivity(frequency)

        start, end = self.__mri.bounding_box()
        return VoxelSpace(data=data, start=tuple(start), end=tuple(end))

    def permitivity(self, frequency: float) -> VoxelSpace:
        csf_position = self.__material_distribution(Material.CSF)
        gm_position = self.__material_distribution(Material.GRAY_MATTER)
        wm_position = self.__material_distribution(Material.WHITE_MATTER)

        csf_model = DielectricModel1.create_model(Material.CSF)
        gm_model = DielectricModel1.create_model(Material.GRAY_MATTER)
        wm_model = DielectricModel1.create_model(Material.WHITE_MATTER)

        default = gm_model.permitivity(frequency)
        data = np.full(self.__mri.data_map().shape, default)
        data[csf_position] = csf_model.permitivity(frequency)
        data[gm_position] = gm_model.permitivity(frequency)
        data[wm_position] = wm_model.permitivity(frequency)

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


@dataclass
class DielectricValues:
    start: tuple
    ent: tuple
    permitivity: np.ndarray
    conductivity: np.ndarray

    def cf_conductivity(self, complex: bool = False):
        type = 'complex' if complex else 'float'
        return ngsolve.VoxelCoefficient(start=self.start,
                                        end=self.end,
                                        values=self.conductivity.astype(type),
                                        linear=False)

    def cf_permitivity(self, complex: bool = False):
        type = 'complex' if complex else 'float'
        return ngsolve.VoxelCoefficient(start=self.start,
                                        end=self.end,
                                        values=self.permitivity.astype(type),
                                        linear=False)
