
from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np


class Signal(ABC):
    """Template for Signals."""

    SPACING_FACTOR = 1e4

    @property
    @abstractmethod
    def frequency(self) -> float:
        pass

    @abstractmethod
    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        pass

    def fft_analysis(self) -> Tuple[np.ndarray]:
        sample_spacing = 1 / (self.frequency * self.SPACING_FACTOR)
        samples = self.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        return complex_values, frequencies


class RectangleSignal(Signal):
    """Creates rectangular signal.

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


class TriangleSignal(Signal):
    """Creates triangular signal.

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


class TrapzoidSignal(Signal):
    """Creates trapzoid signal.

    Parameters
    ----------
    frequency : float
        Frequency [Hz] of the signal.
    pulse_width : float
        Relative pulse width of one period.
    top_width: float
        Relative top width of trapzoid signal.
    counter_pulse_width: float
        Relative width of counter pulse of one period.
    inter_pulse_width: float
        Relative width between pulse and counter pulse of one period.
    """
    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 top_width: float,
                 counter_pulse_width: float,
                 inter_pulse_width: float) -> None:
        self.__frequency = abs(frequency)
        self.__pulse_width = abs(pulse_width)
        self.__top_width = abs(top_width)
        self.__space_width = abs(inter_pulse_width)
        self.__counter_pulse_width = abs(counter_pulse_width)

    @property
    def frequency(self):
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

        spacing = min(abs(sample_spacing), 1) * self.__frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self.__pulse_width * n_samples)
        pulse_top_length = int(self.__top_width * n_samples)
        pulse = self.__create_pulse(pulse_length, pulse_top_length)

        space_length = int(self.__space_width * n_samples)
        counter_length = int(self.__counter_pulse_width * n_samples)
        f = self.__top_width / self.__pulse_width if self.__pulse_width else 0
        counter_top_length = int(f * self.__counter_pulse_width * n_samples)
        counter_value = -pulse_length / counter_length if counter_length else 0
        counter_pulse = counter_value * self.__create_pulse(counter_length,
                                                            counter_top_length)
        counter_start = pulse_length + space_length
        counter_end = counter_start + counter_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        n_counter_samples = n_samples - counter_start
        signal[counter_start:counter_end] = counter_pulse[:n_counter_samples]
        return signal

    @staticmethod
    def __create_pulse(pulse_length: int, pulse_top_length: int) -> np.ndarray:

        if not pulse_top_length:
            return np.array([0] * pulse_length)

        ramp_length = int((pulse_length - pulse_top_length) * 0.5)
        step_size = 1 / (ramp_length + 1)
        ramp = np.arange(step_size, 1, step_size)[:ramp_length]
        return np.concatenate((ramp, [1] * pulse_top_length, np.flip(ramp)))
