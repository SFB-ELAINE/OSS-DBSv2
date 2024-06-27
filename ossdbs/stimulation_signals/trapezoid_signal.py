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

    TODO implement it
    """

    def __init__(
        self,
        frequency: float,
        pulse_width: float,
        top_width: float,
        counter_pulse_width: float,
        inter_pulse_width: float,
        counter_pulse_amplitude: float,
    ) -> None:
        self._top_width = top_width
        super().__init__(
            frequency,
            pulse_width,
            inter_pulse_width,
            counter_pulse_width,
            counter_pulse_amplitude,
        )

    def get_time_domain_signal(self, dt: float, timesteps: int) -> np.ndarray:
        """Build time domain signal."""
        if np.isclose(dt, 0.0):
            return np.array([0])
        if np.isclose(self._top_width, 0.0):
            return np.zeros(timesteps)

        signal = np.zeros(timesteps)
        period = 1.0 / self.frequency
        # use offset for visualization
        offset = round(self._pulse_width / dt)
        while offset < timesteps:
            self._generate_time_domain_signal(
                signal, offset, self._pulse_width, self._top_width, self.amplitude, dt
            )
            if not np.isclose(self._counter_pulse_width, 0.0):
                counter_pulse_start_idx = offset + round(
                    self._pulse_width / dt + self._inter_pulse_width / dt
                )
                self._generate_time_domain_signal(
                    signal,
                    counter_pulse_start_idx,
                    self._counter_pulse_width,
                    self._top_width,
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
        top_width: float,
        amplitude: float,
        dt: float,
    ):
        pulse_duration = round(width / dt)
        top_duration = round(top_width / dt)
        pulse_slope_duration = (pulse_duration - top_duration) // 2

        for i in range(pulse_slope_duration):
            if start_index + i < len(signal):
                signal[start_index + i] = amplitude * (i / pulse_slope_duration)

        for i in range(top_duration):
            if start_index + pulse_slope_duration + i < len(signal):
                signal[start_index + pulse_slope_duration + i] = amplitude

        for i in range(pulse_slope_duration):
            if start_index + pulse_slope_duration + top_duration + i < len(signal):
                signal[
                    start_index + pulse_slope_duration + top_duration + i
                ] = amplitude * (1 - (i / pulse_slope_duration))
