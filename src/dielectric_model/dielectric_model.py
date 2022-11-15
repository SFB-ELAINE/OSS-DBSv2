from abc import ABC, abstractmethod
from src.brainsubstance import Material


class AbstractDielectricModel(ABC):

    @abstractmethod
    def permitivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def conductivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def create_model(cls, material: Material) \
            -> 'AbstractDielectricModel':
        pass
