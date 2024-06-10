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

    def get_fourier_coefficients(self, frequency_list: np.ndarray) -> np.ndarray:
        """Get coefficients of Fourier series."""
        coefficients = self._harmonics_at_freqs(
            frequency_list,
            self.amplitude,
            self.frequency,
            self._pulse_width,
            shift=self._pulse_width,
        )
        if np.greater(self._counter_pulse_width, 0.0):
            coefficients2 = self._harmonics_at_freqs(
                frequency_list,
                -self.counter_amplitude,
                self.frequency,
                self._counter_pulse_width,
                shift=2.0 * self._pulse_width + self._inter_pulse_width,
            )
            coefficients = coefficients + coefficients2
        return coefficients

    def _harmonics_at_freqs(self, frequencies, amp, frequency, tp, shift=None):
        """Get Fourier harmonics."""
        coefficient = amp * tp * frequency * np.sinc(frequencies * tp)
        coefficient = coefficient.astype("complex128")
        if shift is not None:
            # apply time shift at each frequency
            coefficient = coefficient * np.exp(-frequencies * 1j * 2.0 * np.pi * shift)
        return coefficient

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
