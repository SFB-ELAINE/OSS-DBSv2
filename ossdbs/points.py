

from abc import ABC, abstractmethod
import numpy as np
from ossdbs.spectrum_modes.spectrum_point_analysis import TimeResult

from ossdbs.vta_points import VTAPointMatrix


class Points(ABC):

    @abstractmethod
    def coordinates(self) -> np.ndarray:
        pass

    @abstractmethod
    def save(self, time_result: TimeResult, path: str) -> None:
        pass


class VTAPoints(Points):

    def __init__(self, coordinates: np.ndarray) -> None:
        self.__coordinates = coordinates

    def coordinates(self) -> np.ndarray:
        return self.__coordinates

    def save(self, time_result: TimeResult, path: str) -> None:
        time_result.save(path)
