from dataclasses import dataclass
from abc import ABC, abstractmethod
from typing import Tuple
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

    The spectrum is also received from here.

    """

    # TODO check for varying frequencies, pulse widths, non-equidistant samples etc
    SAMPLES_PER_PERIODE = 1e4

    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 inter_pulse_width: float,
                 counter_pulse_width: float = None,
                 ) -> None:
        if np.isclose(frequency, 0):
            raise ValueError("Frequency must be greater than zero.")
        self._frequency = abs(frequency)
        self._pulse_width = abs(pulse_width)
        self._space_width = abs(inter_pulse_width)
        self._counter_pulse_width = abs(counter_pulse_width)

        self._frequencies = None
        self._fourier_coefficients = None
        self._sample_spacing = 1 / (self.frequency * self.SAMPLES_PER_PERIODE)

        # sets values for frequencies and fourier coefficients
        self._fft_spectrum()

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

    @property
    def sample_spacing(self):
        return self._sample_spacing

    @property
    def frequencies(self) -> np.ndarray:
        return self._frequencies

    @property
    def fourier_coefficients(self) -> np.ndarray:
        return self.fourier_coefficients

    @abstractmethod
    def _generate_samples(self, sample_spacing: float) -> np.ndarray:
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
        pass

    def _fft_spectrum(self) -> Tuple[np.ndarray]:
        """Returns the complex values and the frequencies from the FFT of this
        signal.

        Returns
        -------
        tuple of np.ndarray
            First value is the collection of complex values, second value is
            collection of the corresponding frequencies.
        """
        samples = self._generate_samples(self.sample_spacing)
        frequencies = np.fft.rfftfreq(int(self.SAMPLES_PER_PERIODE), self.sample_spacing)
        self._frequencies = frequencies
        self._fourier_coefficients = np.fft.rfft(samples)

    def get_octave_band_spectrum(self,
                                 cutoff_frequency: float = None):
        """TODO document

        """
        n_octaves = int(np.log2(len(self.frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        # TODO check
        # old version
        # octave_frequencies = self.frequency * octave_indices
        octave_frequencies = self.frequencies[octave_indices]
        octave_amplitudes = self.fourier_coefficients[octave_indices]
        # TODO implement coefficients
        if cutoff_frequency is None:
            return octave_frequencies, octave_amplitudes
        return self._truncate_spectrum(cutoff_frequency, octave_frequencies, octave_amplitudes)

    def get_truncated_spectrum(self,
                               cutoff_frequency: float):
        return self._truncate_spectrum(cutoff_frequency, self.frequencies, self.fourier_coefficients)

    def _truncate_spectrum(self,
                           cutoff_frequency: float,
                           frequencies: np.ndarray,
                           coefficients: np.ndarray):
        if cutoff_frequency < frequencies[0]:
            raise ValueError("cutoff_frequency is smaller than lowest frequency in spectrum")
        cutoff_idx = np.where(frequencies < cutoff_frequency)[0][-1]
        return frequencies[:cutoff_idx], coefficients[:cutoff_idx]
