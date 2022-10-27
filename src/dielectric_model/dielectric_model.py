from abc import ABC, abstractmethod
from src.brainsubstance import BrainSubstance


class DielectricModel(ABC):

    @abstractmethod
    def permitivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def conductivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def create_model(cls, material: BrainSubstance) -> 'DielectricModel':
        pass


class DielectricModelCSF(DielectricModel):

    def permitivity(self, frequency: float) -> float:
        return 80

    def conductivity(self, frequency: float) -> float:
        return 1.79

    @classmethod
    def create_model(cls, material: BrainSubstance) -> 'DielectricModel':
        return cls()
