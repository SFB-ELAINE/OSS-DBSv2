

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Tuple
import numpy as np
import h5py
from ossdbs.point_analysis.time_results import TimeResult


class PointModel(ABC):

    @abstractmethod
    def coordinates(self) -> np.ndarray:
        pass

    @abstractmethod
    def save(self, time_result: TimeResult, path: str) -> None:
        pass

    @abstractmethod
    def set_location_names(self, names: np.ndarray) -> None:
        pass
