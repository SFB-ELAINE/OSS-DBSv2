
from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np


class Signal(ABC):
    """Template for Signals."""

    SAMPLES_PER_PERIODE = 1e4

    @property
    @abstractmethod
    def frequency(self) -> float:
        """Return frequency of signal.

        Returns
        -------
        float
        """
        pass

    @abstractmethod
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
        pass

    def fft_analysis(self) -> Tuple[np.ndarray]:
        """Returns the complex values and the frequncies from the FFT of this
        signal.

        Returns
        -------
        tuple of np.ndarray
            First value is the collection of complex values, second value is
            collection of the corresponding frequencies.
        """
        sample_spacing = 1 / (self.frequency * self.SAMPLES_PER_PERIODE)
        samples = self.generate_samples(sample_spacing)
        return np.fft.rfft(samples)

    def fft_frequncies(self) -> np.ndarray:
        sample_spacing = 1 / (self.frequency * self.SAMPLES_PER_PERIODE)
        return np.fft.rfftfreq(int(self.SAMPLES_PER_PERIODE), sample_spacing)
