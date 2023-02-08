
import numpy as np


class RectangleStimulationSignal:

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 counter_pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width
        self.space_width = inter_pulse_width
        self.counter_pulse_width = counter_pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.pulse_width * n_samples)
        pulse = np.array([1] * pulse_length)

        space_length = int(self.space_width * n_samples)

        counter_length = int(self.counter_pulse_width * n_samples)
        counter_value = -pulse_length / counter_length if counter_length else 0
        counter_pulse = np.array([counter_value] * counter_length)
        counter_start = pulse_length + space_length
        counter_end = counter_start + counter_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        n_counter_samples = n_samples - counter_start
        signal[counter_start:counter_end] = counter_pulse[:n_counter_samples]
        return signal


class TriangleStimulationSignal:

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 counter_pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.frequency = frequency
        self.pulse_width = pulse_width
        self.space_width = inter_pulse_width
        self.counter_pulse_width = counter_pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        spacing = sample_spacing * self.frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.pulse_width * n_samples)
        pulse = self.__create_pulse(pulse_length)

        space_length = int(self.space_width * n_samples)

        counter_length = int(self.counter_pulse_width * n_samples)
        counter_value = -pulse_length / counter_length if counter_length else 0
        counter_pulse = counter_value * self.__create_pulse(counter_length)
        counter_start = pulse_length + space_length
        counter_end = counter_start + counter_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        n_counter_samples = n_samples - counter_start
        signal[counter_start:counter_end] = counter_pulse[:n_counter_samples]

        return signal

    @staticmethod
    def __create_pulse(pulse_length: int):
        if not pulse_length:
            return np.array([])

        is_odd = pulse_length % 2
        n_ramp_samples = max(0, int(pulse_length - 2 + is_odd) // 2)
        step_size = 1 / (n_ramp_samples + 1)
        ramp = np.arange(step_size, 1, step_size)[:n_ramp_samples]
        rest_samples = [0] * int(pulse_length - 2 * len(ramp) - 1)
        return np.concatenate((ramp, [1], np.flip(ramp), rest_samples))
