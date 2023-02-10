
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import h5py


@dataclass
class Result:

    points: np.ndarray[float]
    frequency: np.ndarray[float]
    potential: np.ndarray[complex]
    current_density: np.ndarray[complex]
    conductivity: np.ndarray[complex]

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset("points", data=self.points)
            file.create_dataset("frequencies", data=self.frequency)
            file.create_dataset("potential", data=self.potential)
            file.create_dataset("current_density", data=self.current_density)
            file.create_dataset("conductivity", data=self.conductivity)


@dataclass
class ResultImpedance:

    points: np.ndarray[float]
    frequency: np.ndarray[float]
    potential: np.ndarray[complex]
    current_density: np.ndarray[complex]
    conductivity: np.ndarray[complex]
    imdedance: np.ndarray[complex]

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset("points", data=self.points)
            file.create_dataset("frequencies", data=self.frequency)
            file.create_dataset("potential", data=self.potential)
            file.create_dataset("current_density", data=self.current_density)
            file.create_dataset("conductivity", data=self.conductivity)
            file.create_dataset("impedance", data=self.imdedance)


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> None:
        pass
