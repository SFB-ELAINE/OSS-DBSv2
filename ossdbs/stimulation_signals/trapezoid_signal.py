# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from .signal import TimeDomainSignal


class TrapezoidSignal(TimeDomainSignal):
    """Represents trapezoid signal.

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

    def __init__(
        self,
        frequency: float,
        pulse_width: float,
        top_width: float,
        counter_pulse_width: float,
        inter_pulse_width: float,
    ) -> None:
        self._top_width = top_width
        super().__init__(frequency, pulse_width, inter_pulse_width, counter_pulse_width)

    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        """Obtain Fourier coefficients of signal."""
        raise NotImplementedError(
            """Fourier coefficients for
                                     trapezoidal signal not yet implemented."""
        )
