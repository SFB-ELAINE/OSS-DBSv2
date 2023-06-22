
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
        pass

    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        pass
