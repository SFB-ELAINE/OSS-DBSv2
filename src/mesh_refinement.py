
from src.brainsubstance import Material
from src.mesh import Mesh
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
import numpy as np


class MeshRefinement:

    def __init__(self, mesh: Mesh) -> None:
        self.mesh = mesh

    def refine_by_mri(self, mri: MagneticResonanceImage) -> None:
        maximum_size = min(mri.voxel_size())
        csf_voxel = mri.material_distribution(Material.CSF)
        flags = np.logical_and(self.mesh.elements_at_position(csf_voxel),
                               self.mesh.elements_grater_than(maximum_size))
        while np.any(flags) and self.mesh.sobolev_space().ndof < 5e5:
            self.mesh.set_refinement_flag(flags)
            self.mesh.refine()
            csf_voxel = mri.material_distribution(Material.CSF)
            flags = np.logical_and(self.mesh.elements_at_position(csf_voxel),
                                   self.mesh.elements_grater_than(maximum_size)
                                   )

    def refine_by_boundaries(self, boundaries: list) -> None:
        self.mesh.set_boundary_refinement_flags(boundaries)
        self.mesh.refine()
