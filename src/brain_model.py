from src.electrodes import AbstractElectrode
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brain_imaging.diffusion_tensor_imaging import DiffusionTensorImage
from src.mesh import Mesh
import netgen


class BrainModel:

    def __init__(self,
                 mri: MagneticResonanceImage,
                 dti: DiffusionTensorImage = None,
                 electrode: AbstractElectrode = None) -> None:
        self.__mri = mri
        self.__dti = dti
        self.__electrode = electrode

    def generate_mesh(self, order: int = 2, boundary_values: list = None):
        return Mesh(self.generate_geometry,
                    order=order,
                    boundaries=boundary_values)

    def generate_geometry(self):
        if not self.electrode:
            return self.__brain_geometry()
        return self.__brain_geometry - self.__electrode.generate_geometry()

    def __brain_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        x, y, z = self.__mri.bounding_box()[1]
        matrix = [x, 0, 0, 0, y, 0, 0, 0, z]
        trasformator = netgen.occ.gp_GTrsf(mat=matrix)
        return trasformator(netgen.occ.Sphere((1, 1, 1), 1))
