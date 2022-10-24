from dataclasses import dataclass
from src.Nifti1Image import Nifti1Image
from src.brainsubstance import BrainSubstance


@dataclass
class MaterialMap:
    csf_index: int = 1
    wm_index: int = 2
    gm_index: int = 3

    def material(self, index: int) -> BrainSubstance:
        if self.csf_index == index:
            return BrainSubstance.CEREBROSPINAL_FLUID

        if self.wm_index == index:
            return BrainSubstance.WHITE_MATTER

        if self.gm_index == index:
            return BrainSubstance.GRAY_MATTER

        return BrainSubstance.UNKNOWN


class MagneticResonanceImage(Nifti1Image):

    def __init__(self,
                 file_path: str,
                 material_map: MaterialMap = MaterialMap()) -> None:
        super().__init__(file_path)
        self.__map = material_map

    # def material_at(self, positions: list) -> list[BrainSubstance]:

    #     return [self.__map.material(value)
    #             for value in self.values_at(positions)]
