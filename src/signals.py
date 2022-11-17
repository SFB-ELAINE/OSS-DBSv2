import numpy as np


class RectangleSignal:

    def __init__(self, frequency: float, pulse_width: float) -> None:
        self.__frequency = frequency
        self.__pulse_width = pulse_width

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        samples = np.zeros(n_samples)
        samples[:int(n_samples * self.__pulse_width)] = 1.0
        return samples


class TriangleSignal:

    def __init__(self, frequency: float, pulse_width: float) -> None:
        self.__frequency = frequency
        self.__pulse_width = pulse_width

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        pulse = self.__create_pulse(n_samples)
        samples = np.zeros(n_samples)
        samples[:len(pulse)] = pulse
        return samples

    def __create_pulse(self, n_samples):
        pulse_samples = n_samples * self.__pulse_width
        is_odd = pulse_samples % 2
        n_ramp_samples = int(pulse_samples - 2 + int(is_odd)) // 2
        step_size = 1 / (n_ramp_samples + 1)
        right_ramp = np.arange(1, 0, -step_size)
        left_ramp = np.flip(right_ramp)
        pulse = np.concatenate((left_ramp, right_ramp[1:]))
        return pulse


class TrapzoidSignal:

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 top_width: float) -> None:
        self.__frequency = frequency
        self.__pulse_width = pulse_width
        self.__top_width = top_width

    def generate_samples(self, sample_spacing: float):
        n_samples = int(1 / (sample_spacing * self.__frequency))
        pulse = self.__create_pulse(n_samples)
        samples = np.zeros(n_samples)
        samples[:len(pulse)] = pulse
        return samples

    def __create_pulse(self, n_samples):
        ramp_length = (self.__pulse_width - self.__top_width) * 0.5
        step_size = 1 / int(n_samples * ramp_length + 1)
        right_ramp = np.arange(1, 0, -step_size)[1:]
        left_ramp = np.flip(right_ramp)
        top_samples = int(self.__pulse_width * n_samples) - 2 * len(right_ramp)
        pulse = np.concatenate((left_ramp, [1] * top_samples, right_ramp))
        return pulse
