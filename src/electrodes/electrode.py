from abc import ABC, abstractmethod


class Electrode(ABC):

    @abstractmethod
    def translate(self, vector: tuple) -> 'Electrode':
        pass

    @abstractmethod
    def rotate(self, angle: float) -> 'Electrode':
        pass

    @abstractmethod
    def tilt(self, direction: tuple) -> 'Electrode':
        pass
