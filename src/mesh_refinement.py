
from src.brainsubstance import Material
from src.mesh import Mesh
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
import numpy as np


class MeshRefinement:

    def __init__(self, mesh: Mesh) -> None:
        self.mesh = mesh

    def refine_by_mri(self, mri: MagneticResonanceImage) -> Mesh:
        flags = self.__elements_to_refine(mri)

        while np.any(flags):
            self.mesh.set_refinement_flag(flags)
            self.mesh.refine()
            flags = self.__elements_to_refine(mri)

    def __elements_to_refine(self, mri):
        maximum_size = min(mri.voxel_size())
        csf_voxel = mri.material_distribution(Material.CSF)
        return np.logical_and(self.mesh.elements_at_position(csf_voxel),
                              self.mesh.elements_grater_than(maximum_size))
