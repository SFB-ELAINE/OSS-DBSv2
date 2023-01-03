from abc import ABC, abstractmethod
import numpy as np


class Signal(ABC):

    @abstractmethod
    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        pass


class RectangleSignal(Signal):

    def __init__(self, frequency: float, pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1
        pulse = self.__create_pulse(n_samples)
        padding = np.zeros(n_samples - len(pulse))
        return np.concatenate((pulse[:n_samples], padding))

    def __create_pulse(self, n_samples: int) -> np.ndarray:
        pulse_samples = int(n_samples * min(self.pulse_width, 1))
        return np.array([1] * pulse_samples)


class TriangleSignal(Signal):

    def __init__(self, frequency: float, pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1
        pulse = self.__create_pulse(n_samples)
        padding = np.zeros(n_samples - len(pulse))
        return np.concatenate((pulse, padding))

    def __create_pulse(self, n_samples: int) -> np.ndarray:
        pulse_samples = int(n_samples * min(self.pulse_width, 1))
        is_odd = pulse_samples % 2
        n_ramp_samples = max(0, int(pulse_samples - 2 + is_odd) // 2)
        step_size = 1 / (n_ramp_samples + 1)
        ramp = np.arange(step_size, 1, step_size)[:n_ramp_samples]
        peak_sample = pulse_samples > 0
        return np.concatenate((ramp, [1] * peak_sample, np.flip(ramp)))


class TrapzoidSignal(Signal):

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 top_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width
        self.top_width = top_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1
        pulse = self.__create_pulse(n_samples)
        padding = np.zeros(n_samples - len(pulse))
        return np.concatenate((pulse, padding))

    def __create_pulse(self, n_samples: int) -> np.ndarray:
        pulse_width = min(self.pulse_width, 1)
        ramp_length = (pulse_width - self.top_width) * 0.5
        n_ramp_samples = max(0, int(n_samples * ramp_length))
        step_size = 1 / int(n_ramp_samples + 1)
        ramp = np.arange(step_size, 1, step_size)[:n_ramp_samples]
        top_samples = int(pulse_width * n_samples) - 2 * len(ramp)
        return np.concatenate((ramp, [1] * top_samples, np.flip(ramp)))
