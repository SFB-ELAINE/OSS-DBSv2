
from ossdbs.stimmulation_signals.signal import Signal
import numpy as np


class RectangleSignal(Signal):
    """Represents a rectangular signal.

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
        pulse = np.array([1] * pulse_length)

        space_length = int(self.__space_width * n_samples)

        counter_length = int(self.__counter_pulse_width * n_samples)
        counter_value = -pulse_length / counter_length if counter_length else 0
        counter_pulse = np.array([counter_value] * counter_length)
        counter_start = pulse_length + space_length
        counter_end = counter_start + counter_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        n_counter_samples = n_samples - counter_start
        signal[counter_start:counter_end] = counter_pulse[:n_counter_samples]
        return signal
