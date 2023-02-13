
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import h5py


@dataclass
class TimeResult:
    points: np.ndarray[float]
    time_steps: np.ndarray[float]
    potential: np.ndarray[float]
    current_density: np.ndarray[float]


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

    def time_result(self):

        n_time_steps = (len(self.frequency) - 1) * 2
        potential_t = np.zeros((len(self.points), n_time_steps))
        current_density_t = np.zeros((len(self.points), n_time_steps, 3))

        for start in range(0, len(self.points), 1000):
            end = start + 1000
            potential = self.potential[start:end]
            potential_t[start:end] = np.fft.irfft(potential, axis=1)
            current_density = self.current_density[start:end]
            current_density_t[start:end] = np.fft.irfft(current_density,
                                                        axis=1)

        frequency = self.frequency[1]
        n_steps = self.potential.shape[1]
        time = np.arange(n_steps) * 1 / (frequency * n_steps)
        return TimeResult(points=self.points,
                          potential=potential_t,
                          current_density=current_density_t,
                          time_steps=time)


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> None:
        pass
