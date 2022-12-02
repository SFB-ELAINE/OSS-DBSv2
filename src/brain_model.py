
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brainsubstance import Material
from src.dielectric_model.dielectric_model_var1 import ModelCreator
from src.electrodes.abstract_electrode import AbstractElectrode
from src.geometry import Geometry
from src.mesh import Mesh
import netgen
import numpy as np

from src.voxels import Voxels


class BrainModel:

    __MATERIALS = [Material.CSF, Material.GRAY_MATTER, Material.WHITE_MATTER]

    def __init__(self, mri: MagneticResonanceImage) -> None:
        self.__mri = mri
        self.__electrodes = []

    def complex_conductivity(self, frequency: float) -> Voxels:
        omega = 2 * np.pi * frequency
        default = ModelCreator.create(Material.GRAY_MATTER).conductivity(omega)
        data = np.full(self.__mri.xyz_shape(), default)

        for material in self.__MATERIALS:
            position = self.__mri.material_distribution(material=material)
            conductivity = ModelCreator.create(material).conductivity(omega)
            data[position.data] = conductivity

        start, end = self.__mri.bounding_box()
        return Voxels(data=data, start=start, end=end)

    def add_electrode(self, electrode: AbstractElectrode):
        self.__electrodes.append(electrode)

    def generate_mesh(self, order: int = 2):
        return Mesh(self.generate_geometry(), order=order)

    def generate_geometry(self):
        brain = self.__create_ellipsoid()

        geometry = brain
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()

        return Geometry(geometry=geometry)

    def __create_ellipsoid(self):
        start, end = self.__mri.bounding_box()
        x, y, z = (np.array(end) - np.array(start)) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = trasformator(sphere).Move(tuple(start))
        ellipsoid.bc('Brain')
        return ellipsoid
