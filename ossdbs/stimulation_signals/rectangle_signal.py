# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from .signal import TimeDomainSignal


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

    def get_time_domain_signal(self, dt: float, timesteps: int) -> np.ndarray:
        """Build time domain signal."""
        signal = np.zeros(timesteps)
        period = 1.0 / self.frequency
        # use offset for visualization
        offset = int(self._pulse_width / dt)
        while offset < timesteps:
            pulse_idx = offset + int(self._pulse_width / dt)
            signal[offset:pulse_idx] = self.amplitude
            if not np.isclose(self._counter_pulse_width, 0.0):
                counter_pulse_start_idx = offset + int(
                    self._pulse_width / dt + self._inter_pulse_width / dt
                )
                counter_pulse_end_idx = counter_pulse_start_idx + int(
                    self._counter_pulse_width / dt
                )
                signal[
                    counter_pulse_start_idx:counter_pulse_end_idx
                ] = -self.counter_amplitude
            offset += int(period / dt)
        return signal
