import numpy as np

from .signal import TimeDomainSignal


class TriangleSignal(TimeDomainSignal):
    """Represents triangular signal."""

    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        """Get coefficients of Fourier series."""
        pass
