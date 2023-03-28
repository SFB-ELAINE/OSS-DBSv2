
from ossdbs.stimmulation_signals.signal import Signal
import numpy as np


class TriangleSignal(Signal):
    """Represents triangular signal.

    Parameters
    ----------
    frequency : float
        Frequency [Hz] of the signal.
    pulse_width : float
        Relative pulse width of one period.
    counter_pulse_width: float
        Relative width of counter pulse of one period.
    inter_pulse_width: float
        Relative width between pulse and counter pulse of one period.
    """
    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 counter_pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.__frequency = abs(frequency)
        self.__pulse_width = abs(pulse_width)
        self.__space_width = abs(inter_pulse_width)
        self.__counter_pulse_width = abs(counter_pulse_width)

    @property
    def frequency(self) -> float:
        """Return frequency of signal.

        Returns
        -------
        float
        """
        return self.__frequency

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        """Generate samples which follow the signal form.

        Parameters
        ----------
        sample_spacing : float
            Timestep [s] between two samples.

        Returns
        -------
        np.ndarray
            Samples for one period.
        """

        spacing = sample_spacing * self.__frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.__pulse_width * n_samples)
        pulse = self.__create_pulse(pulse_length)

        space_length = int(self.__space_width * n_samples)

        counter_length = int(self.__counter_pulse_width * n_samples)
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
    def __create_pulse(pulse_length: int) -> np.ndarray:
        is_odd = pulse_length % 2
        n_ramp_samples = max(0, int(pulse_length - 2 + is_odd) // 2)
        step_size = 1 / (n_ramp_samples + 1)
        ramp = np.arange(step_size, 1, step_size)[:n_ramp_samples]
        rest_samples = [0] * int(pulse_length - 2 * len(ramp) - 1)
        return np.concatenate((ramp, [1], np.flip(ramp), rest_samples))
