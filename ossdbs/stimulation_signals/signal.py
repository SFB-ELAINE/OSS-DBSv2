# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq

from .utilities import adjust_cutoff_frequency, retrieve_time_domain_signal_from_fft


@dataclass
class FrequencyDomainSignal:
    """Store information for freqency domain signal."""

    frequencies: np.ndarray
    amplitudes: np.ndarray
    current_controlled: bool
    base_frequency: float
    cutoff_frequency: float
    signal_length: int
    octave_band_approximation: bool = False


class TimeDomainSignal(ABC):
    """Template for Signals.

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

    Notes
    -----
    TODO document and clarify how to use amplitude
    The spectrum is also received from here.

    """

    def __init__(
        self,
        frequency: float,
        pulse_width: float,
        inter_pulse_width: float,
        counter_pulse_width: Optional[float] = None,
        counter_pulse_amplitude: Optional[float] = 1.0,
    ) -> None:
        if np.isclose(frequency, 0):
            raise ValueError("Frequency must be greater than zero.")
        self._frequency = frequency
        self._pulse_width = pulse_width
        self._inter_pulse_width = inter_pulse_width
        self._counter_pulse_width = counter_pulse_width
        # the values here are relative amplitudes
        # e.g., if the signal is 1V and the counter_pulse_amplitude is 0.5
        # its amplitude will be 0.5V (same for 1mA)
        self._amplitude = 1.0
        self._counter_amplitude = counter_pulse_amplitude

    @property
    def amplitude(self) -> float:
        """Return signal amplitude."""
        return self._amplitude

    @amplitude.setter
    def amplitude(self, value) -> None:
        """Set amplitude value."""
        self._amplitude = value

    @property
    def counter_amplitude(self) -> float:
        """Get amplitude of counterpulse."""
        return self._counter_amplitude

    @counter_amplitude.setter
    def counter_amplitude(self, value) -> None:
        self._counter_amplitude = value

    @property
    def frequency(self) -> float:
        """Return frequency of signal.

        Returns
        -------
        float
        """
        return self._frequency

    @frequency.setter
    def frequency(self, value):
        """Set frequency of signal."""
        self._frequency = value

    @abstractmethod
    def get_fourier_coefficients(frequencies: float) -> np.ndarray:
        """Obtain Fourier coefficients of signal."""
        pass

    def get_frequencies_and_fourier_coefficients(
        self, cutoff_frequency: float
    ) -> np.ndarray:
        """Get frequencies and Fourier coefficients up to cutoff frequency.

        Parameters
        ----------
        cutoff_frequency: float
            Highest considered frequency.
        """
        max_harmonic = int(cutoff_frequency / self.frequency)
        harmonics = np.arange(0, max_harmonic + 1)
        frequencies = harmonics * self.frequency
        coefficients = self.get_fourier_coefficients(frequencies)
        return frequencies, coefficients

    def get_fft_spectrum(self, cutoff_frequency: float) -> np.ndarray:
        """FFT spectrum of time-domain signal.

        Parameters
        ----------
        cutoff_frequency: float
            Highest considered frequency.

        """
        # we double the cutoff_frequency to actually sample until there
        cutoff_frequency = adjust_cutoff_frequency(
            2.0 * cutoff_frequency, self.frequency
        )
        dt = 1.0 / cutoff_frequency
        # required length for frequency
        timesteps = int(cutoff_frequency / self.frequency)

        time_domain_signal = self.get_time_domain_signal(dt, timesteps)
        return fftfreq(len(time_domain_signal), d=dt), fft(time_domain_signal)

    def retrieve_time_domain_signal(
        self, fft_signal: np.ndarray, cutoff_frequency: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute time-domain signal by FFT."""
        return retrieve_time_domain_signal_from_fft(
            fft_signal, cutoff_frequency, self.frequency
        )

    @abstractmethod
    def get_time_domain_signal(self, dt: float, timesteps: int) -> np.ndarray:
        """Time-domain signal for given timestep."""
        pass

    def plot_time_domain_signal(self, cutoff_frequency, output_path, show=False):
        """Plot signal and export to PDF."""
        cutoff_frequency = adjust_cutoff_frequency(
            2.0 * cutoff_frequency, self.frequency
        )
        dt = 1.0 / cutoff_frequency
        # required length for frequency
        timesteps = int(cutoff_frequency / self.frequency)
        time_domain_signal = self.get_time_domain_signal(dt, timesteps)
        plt.plot(dt * np.arange(0, timesteps), time_domain_signal)
        plt.xlabel("Time / s")
        plt.ylabel("Signal / arb. u.")
        plt.savefig(os.path.join(output_path, "time_domain_signal.pdf"))
        if show:
            plt.show()
        else:
            plt.close()
