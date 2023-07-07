from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np


@dataclass
class FrequencyDomainSignal:
    frequencies: np.ndarray
    amplitudes: np.ndarray
    current_controlled: bool


class TimeDomainSignal(ABC):
    """Template for Signals.

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

    Notes
    -----

    TODO document and clarify how to use amplitude
    The spectrum is also received from here.

    """

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 inter_pulse_width: float,
                 counter_pulse_width: float = None,
                 ) -> None:
        if np.isclose(frequency, 0):
            raise ValueError("Frequency must be greater than zero.")
        self._frequency = frequency
        self._pulse_width = pulse_width
        self._inter_pulse_width = inter_pulse_width
        self._counter_pulse_width = counter_pulse_width
        self._amplitude = 1.0
        self._counter_amplitude = 1.0

    @property
    def amplitude(self) -> float:
        return self._amplitude

    @amplitude.setter
    def amplitude(self, value) -> None:
        self._amplitude = value

    @property
    def counter_amplitude(self) -> float:
        return self._counter_amplitude

    @counter_amplitude.setter
    def counter_amplitude(self, value) -> None:
        self._counter_amplitude = value

    @property
    def frequency(self) -> float:
        """Return frequency of signal.

        Returns
        -------
        float
        """
        return self._frequency

    @frequency.setter
    def frequency(self, value):
        self._frequency = value

    @abstractmethod
    def get_frequencies_and_fourier_coefficients(self,
                                                 cutoff_frequency: float) -> np.ndarray:
        pass

    @abstractmethod
    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        pass

    def get_octave_band_spectrum(self,
                                 cutoff_frequency: float) -> np.ndarray:
        """TODO document

        """
        frequencies, fourier_coefficients = self.get_frequencies_and_fourier_coefficients(cutoff_frequency)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        # TODO check
        # old version
        octave_frequencies = frequencies[octave_indices]
        octave_amplitudes = fourier_coefficients[octave_indices]
        return octave_frequencies, octave_amplitudes
