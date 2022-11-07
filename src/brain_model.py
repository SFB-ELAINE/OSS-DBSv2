from src.dielectric_model import DielectricModel1
from src.brainsubstance import Material
from src.electrodes import AbstractElectrode
from src.brain_imaging import MagneticResonanceImage
from src.brain_imaging import DiffusionTensorImage
from src.mesh import Mesh
import netgen
import numpy as np


class BrainModel:

    def __init__(self,
                 mri: MagneticResonanceImage,
                 dti: DiffusionTensorImage = None,
                 electrode: AbstractElectrode = None) -> None:
        self.__mri = mri
        self.__dti = dti
        self.__electrode = electrode

    def bounding_box(self) -> tuple:
        return self.__mri.bounding_box()

    def material_distribution(self, material: Material) -> np.ndrray:
        return self.__mri.data_map == material

    def conductivity(self, frequency: float) -> None:
        csf_position = self.material_distribution(Material.CSF)
        gm_position = self.material_distribution(Material.GRAY_MATTER)
        wm_positiom = self.material_distribution(Material.WHITE_MATTER)
        unknown_position = self.material_distribution(Material.UNKNOWN)

        csf_model = DielectricModel1.create_model(Material.CSF)
        gm_model = DielectricModel1.create_model(Material.GRAY_MATTER)
        wm_model = DielectricModel1.create_model(Material.WHITE_MATTER)

        conductivity = np.zeros(self.__mri.data_map.shape)
        conductivity[csf_position] = csf_model.conductivity(frequency)
        conductivity[gm_position] = gm_model.conductivity(frequency)
        conductivity[wm_positiom] = wm_model.conductivity(frequency)
        conductivity[unknown_position] = gm_model.conductivity(frequency)
        return conductivity

    def generate_mesh(self, order: int = 2, boundary_values: list = None):
        return Mesh(self.generate_geometry(self.__electrode),
                    order=order,
                    boundaries=boundary_values)

    def generate_geometry(self, electrode: AbstractElectrode = None):
        x, y, z = self.__mri.bounding_box()[1]
        brain = BrainGeometry(x=x, y=y, z=z)
        if not electrode:
            return brain
        return brain - electrode.generate_geometry()


class BrainGeometry:

    def __init__(self, x: int, y: int, z: int) -> None:
        self.__x = x
        self.__y = y
        self.__z = z

    def create(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        matrix = [self.__x, 0, 0, 0, self.__y, 0, 0, 0, self.__z]
        trasformator = netgen.occ.gp_GTrsf(mat=matrix)
        return trasformator(netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1))
