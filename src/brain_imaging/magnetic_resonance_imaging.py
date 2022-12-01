from src.brain_imaging.Nifti1Image import Nifti1Image
from src.brainsubstance import Material
from src.voxels import Voxels


class MagneticResonanceImage(Nifti1Image):

    def material_distribution(self, material: Material) -> Voxels:
        start, end = self.bounding_box()
        is_material = self._image.get_fdata() == material
        return Voxels(data=is_material, start=tuple(start), end=tuple(end))
