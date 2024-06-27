# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from .signal import TimeDomainSignal


class TriangleSignal(TimeDomainSignal):
    """Represents a triangular signal.

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
        if np.isclose(dt, 0.0):
            return np.array([0])
        if np.isclose(self._pulse_width, 0.0):
            return np.zeros(timesteps)

        signal = np.zeros(timesteps)
        period = 1.0 / self.frequency
        # use offset for visualization
        offset = round(self._pulse_width / dt)

        while offset < timesteps:
            self._generate_time_domain_signal(
                signal, offset, self._pulse_width, self.amplitude, dt
            )
            if not np.isclose(self._counter_pulse_width, 0.0):
                counter_pulse_start_index = offset + round(
                    self._pulse_width / dt + self._inter_pulse_width / dt
                )
                self._generate_time_domain_signal(
                    signal,
                    counter_pulse_start_index,
                    self._counter_pulse_width,
                    -self._counter_amplitude,
                    dt,
                )

            offset += round(period / dt)
        return signal

    def _generate_time_domain_signal(
        self,
        signal: np.array,
        start_index: int,
        width: float,
        amplitude: float,
        dt: float,
    ):
        pulse_duration = round(width / dt)
        for i in range(pulse_duration):
            if start_index + i < len(signal):
                if i <= pulse_duration / 2:
                    signal[start_index + i] = amplitude * (i / pulse_duration) * 2
                else:
                    signal[start_index + i] = amplitude * (1 - (i / pulse_duration)) * 2
