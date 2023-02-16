
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import h5py


@dataclass
class TimeResult:
    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset('TimeSteps[s]', data=self.time_steps)
            file.create_dataset('Points[mm]', data=self.points)
            file.create_dataset('Potential[V]', data=self.potential)
            file.create_dataset('Current_density[A/m2]',
                                data=self.current_density)

    def save_by_categories(self, path: str, categories: list) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset('TimeSteps', data=self.time_steps)
            start = 0
            for category in categories:
                name, n_points = category
                end = start + n_points
                h5_group = file.create_group(name)
                points = self.points[start:end]
                h5_group.create_dataset('Points', data=points)
                potential = self.potential[start:end]
                h5_group.create_dataset('Potential', data=potential)
                current_density = self.current_density[start:end]
                h5_group.create_dataset('CurrentDensity', data=current_density)


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


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> None:
        pass
