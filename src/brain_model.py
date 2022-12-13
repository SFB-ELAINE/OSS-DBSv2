
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brainsubstance import Material
from src.dielectric_model.dielectric_model_var1 import ModelCreator
from src.electrodes.abstract_electrode import AbstractElectrode
from src.geometry import Geometry
from src.mesh import Mesh
import netgen
import numpy as np
from src.mesh_refinement import MeshRefinement


class BrainModel:

    def __init__(self, bounding_box) -> None:
        self.__bounding_box = bounding_box
        self.__electrodes = []

    def add_electrodes(self, electrodes: list[AbstractElectrode]):
        for electrode in electrodes:
            self.__electrodes.append(electrode)

    def generate_mesh(self, order: int = 2):
        mesh = Mesh(self.__generate_geometry(), order=order)
        return mesh

    def __generate_geometry(self):
        geometry = self.__create_ellipsoid()
        for electrode in self.__electrodes:
            geometry = geometry - electrode.generate_geometry()
        return Geometry(geometry=geometry)

    def __create_ellipsoid(self):
        start, end = self.__bounding_box
        x, y, z = (np.array(end) - np.array(start)) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = trasformator(sphere).Move(start)
        ellipsoid.bc('Brain')
        return ellipsoid
