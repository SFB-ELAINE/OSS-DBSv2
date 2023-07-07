
from .signal import TimeDomainSignal
import numpy as np


class RectangleSignal(TimeDomainSignal):
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

    def get_frequencies_and_fourier_coefficients(self,
                                                 cutoff_frequency: float) -> np.ndarray:
        max_harmonic = int(cutoff_frequency / self.frequency)
        harmonics = np.arange(0, max_harmonic + 1)
        frequencies = harmonics * self.frequency
        coefficients = self.get_fourier_coefficients(frequencies)
        return frequencies, coefficients

        pass

    def get_fourier_coefficients(self, frequency_list: np.ndarray) -> np.ndarray:
        coefficients = self._harmonics_at_freqs(frequency_list, self.amplitude, self.frequency, self._pulse_width)
        if np.greater(self._counter_pulse_width, 0.0):
            coefficients2 = self._harmonics_at_freqs(frequency_list, -self.counter_amplitude, self.frequency, self._counter_pulse_width, shift=self._pulse_width)
            coefficients = coefficients + coefficients2
        return coefficients

    def _harmonics_at_freqs(self, harmonics, amp, frequency, tp, shift=None):
        coefficient = amp * tp * frequency * np.sinc(harmonics * tp * frequency)
        coefficient = coefficient.astype("complex128")
        if shift is not None:
            # apply time shift at each frequency
            coefficient = coefficient * np.exp(-harmonics * 1j * 2. * np.pi * frequency * shift)
        return coefficient
