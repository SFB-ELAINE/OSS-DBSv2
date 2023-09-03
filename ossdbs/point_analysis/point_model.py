
from abc import ABC, abstractmethod
import numpy as np
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.fem import Mesh


class PointModel(ABC):

    @property
    def coordinates(self) -> np.ndarray:
        return self._coordinates

    @abstractmethod
    def _initialize_coordinates(self) -> np.ndarray:
        pass

    @abstractmethod
    def save(self, time_result: TimeResult, path: str) -> None:
        pass

    @abstractmethod
    def set_location_names(self, names: np.ndarray) -> None:
        pass

    def points_in_mesh(self, mesh: Mesh):
        mask = mesh.not_included(self.coordinates)
        return np.ma.masked_array(self.coordinates, np.column_stack((mask, mask, mask)))
