from .signal import TimeDomainSignal
import numpy as np


class TriangleSignal(TimeDomainSignal):
    """Represents triangular signal.
    """
    def get_frequencies_and_fourier_coefficients(self,
                                                 cutoff_frequency: float) -> np.ndarray:
        pass

    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        pass
