import numpy as np


class RectangleSignal:

    def __init__(self, frequency: float, pulswidth: float) -> None:
        self.__frequency = frequency
        self.__pulswidth = pulswidth

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        samples = np.zeros(n_samples)
        samples[:int(n_samples * self.__pulswidth)] = 1.0
        return samples


class TriangleSignal:

    def __init__(self, frequency: float, pulswidth: float) -> None:
        self.__frequency = frequency
        self.__puls_width = pulswidth

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        n_ramp_samples = int(n_samples * self.__puls_width + 1) // 2
        step_size = 1 / n_ramp_samples
        right_ramp = np.arange(1, 0, -step_size)
        left_ramp = np.flip(right_ramp)
        rest_signal = np.zeros(n_samples - 2 * n_ramp_samples + 1)
        samples = np.concatenate((left_ramp[:-1], right_ramp, rest_signal))
        return samples


class TrapzoidSignal:

    def __init__(self,
                 frequency: float,
                 pulswidth: float,
                 top_width: float) -> None:
        self.__frequency = frequency
        self.__puls_width = pulswidth
        self.__top_width = top_width

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        ramp_length = (self.__puls_width - self.__top_width) * 0.5
        n_ramp_samples = int(n_samples * ramp_length) + 1
        top_samples = int(self.__top_width * n_samples) - 2
        step_size = 1 / n_ramp_samples
        right_ramp = np.arange(1, 0, -step_size)
        left_ramp = np.flip(right_ramp)
        rest_signal = np.zeros(int(n_samples * (1 - self.__puls_width)))
        trapezoid = np.concatenate((left_ramp, [1] * top_samples, right_ramp))
        samples = np.concatenate((trapezoid, rest_signal))
        return samples
