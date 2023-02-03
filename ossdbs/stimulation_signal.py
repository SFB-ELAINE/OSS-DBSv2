
import numpy as np


class RectangleStimulationSignal:

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width
        self.space_width = inter_pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.pulse_width * n_samples)
        pulse = np.array([1] * pulse_length)

        space_length = int(self.space_width * n_samples)

        counter_pulse = np.array([-0.1] * pulse_length * 10)
        counter_start = pulse_length + space_length
        counter_end = counter_start + 10 * pulse_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        signal[counter_start:counter_end] = counter_pulse
        return signal


class TriangleStimulationSignal:

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width
        self.space_width = inter_pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.pulse_width * n_samples)
        pulse = self.__create_pulse(pulse_length)

        space_length = int(self.space_width * n_samples)

        counter_pulse = -0.1 * self.__create_pulse(pulse_length * 10)
        counter_start = pulse_length + space_length
        counter_end = counter_start + 10 * pulse_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        signal[counter_start:counter_end] = counter_pulse

        return signal

    @staticmethod
    def __create_pulse(pulse_length: int):
        if not pulse_length:
            return np.array([])

        is_odd = pulse_length % 2
        n_ramp_samples = max(0, int(pulse_length - 2 + is_odd) // 2)
        step_size = 1 / (n_ramp_samples + 1)
        ramp = np.arange(step_size, 1, step_size)[:n_ramp_samples]
        return np.concatenate((ramp, 1, np.flip(ramp)))
