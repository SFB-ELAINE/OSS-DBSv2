from dataclasses import dataclass
from enum import IntEnum, unique
from src.Nifti1Image import Nifti1Image
import numpy as np

@unique
class BrainSubstance(IntEnum):
    UNKNOWN = 0
    CEREBROSPINAL_FLUID = 1
    WHITE_MATTER = 2
    GRAY_MATTER = 3


@dataclass    
class MaterialMap:
    csf_index: int = 1
    wm_index: int = 2
    gm_index: int = 3

    def material(self, index: int):
        if self.csf_index == index:
            return BrainSubstance.CEREBROSPINAL_FLUID
        
        if self.wm_index == index:
            return BrainSubstance.WHITE_MATTER

        if self.gm_index == index:
            return BrainSubstance.GRAY_MATTER

        return BrainSubstance.UNKNOWN
    

class MagneticResonanceImage(Nifti1Image):

    def __init__(self, file_path, material_map: MaterialMap = None) -> None:
        super().__init__(file_path)
        self.__map = material_map
        if not material_map:
            self.__map = MaterialMap()

    def material_at(self, positions: np.ndarray) -> list[str]:
        return [self.__map.material(value) 
                for value in self.values_at(positions)]
