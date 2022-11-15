from abc import ABC, abstractmethod
from ..brainsubstance import Material


class AbstractDielectricModel(ABC):

    @abstractmethod
    def relative_permitivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def conductivity(self, frequency: float) -> float:
        pass

    @abstractmethod
    def create_model(cls, material: Material) \
            -> 'AbstractDielectricModel':
        pass
