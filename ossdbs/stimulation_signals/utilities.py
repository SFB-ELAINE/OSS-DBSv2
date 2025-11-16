# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from typing import Union

import numpy as np
from scipy.fft import ifft

_logger = logging.getLogger(__name__)


def adjust_cutoff_frequency(cutoff_frequency, frequency):
    """Function to make cutoff frequency multiple of stimulation frequency."""
    return cutoff_frequency - cutoff_frequency % frequency


def get_timesteps(
    cutoff_frequency: float, base_frequency: float, n_frequencies: int
) -> np.ndarray:
    """Return list with timesteps."""
    cutoff_frequency = adjust_cutoff_frequency(2.0 * cutoff_frequency, base_frequency)
    dt = 1.0 / cutoff_frequency
    timesteps = dt * np.arange(n_frequencies)
    return timesteps


def retrieve_time_domain_signal_from_fft(
    fft_signal: np.ndarray, cutoff_frequency: float, base_frequency: float
) -> tuple[np.ndarray, np.ndarray]:
    """Compute time-domain signal via fft."""
    # double the cutoff_frequency to actually sample until there
    signal = ifft(fft_signal).real
    timesteps = get_timesteps(cutoff_frequency, base_frequency, len(signal))
    return timesteps, signal


def reconstruct_time_signals(
    freq_domain_signal: np.ndarray, signal_length: int
) -> np.ndarray:
    """Compute time signals from frequency-domain data.

    Parameters
    ----------
    freq_domain_signal: np.ndarray
      Frequency-domain signal to be transformed
    signal_length: int
      Length of initial time-domain signal

    Notes
    -----
    According to the signal length, the highest frequency is
    a negative or positive frequency (as it comes from the FFT).
    """
    # For even signals, highest frequency is not in positive frequencies
    positive_freqs_part = freq_domain_signal
    if signal_length % 2 == 0:
        positive_freqs_part = freq_domain_signal[:-1]

    # Append the reverted signal without the DC frequency
    tmp_freq_domain = np.append(
        positive_freqs_part,
        np.conjugate(np.flip(freq_domain_signal[1:], axis=0)),
        axis=0,
    )
    # run ifft with maximum possible amount of workers
    result_in_time = ifft(tmp_freq_domain, axis=0, workers=-1)
    return result_in_time.real


def get_octave_band_indices(frequencies: np.ndarray) -> np.ndarray:
    """Return indices of octave band frequencies."""
    n_octaves = int(np.log2(len(frequencies) - 1)) + 1
    octave_indices = 2 ** np.arange(0, n_octaves)
    return octave_indices


def get_minimum_octave_band_index(freq_idx: int) -> int:
    """Get index of lowest frequency in octave band."""
    return int(np.round(freq_idx / np.sqrt(2)))


def get_maximum_octave_band_index(freq_idx: int) -> int:
    """Get index of highest frequency in octave band."""
    return int(np.round(freq_idx * np.sqrt(2)))


def get_indices_in_octave_band(
    freq_idx: int, frequency_indices: list, cutoff_frequency_index: int
) -> Union[list, np.ndarray]:
    """Get indices of frequencies in octave band.

    Notes
    -----
    We start evaluating from the bottom. I.e., it is checked if there is
    an overlap with frequencies from the octave band below (already computed).
    The minimum frequencies are increased until there is no overlap with the
    previous band.
    """
    min_freq = get_minimum_octave_band_index(freq_idx)
    max_freq = get_maximum_octave_band_index(freq_idx)
    list_index = np.argwhere(freq_idx == frequency_indices)
    if list_index.shape != (1, 1):
        raise ValueError("Wrong frequencies for band evaluation supplied")
    list_index = list_index[0][0]
    if freq_idx > 0:
        max_of_prev_band = get_maximum_octave_band_index(
            frequency_indices[list_index - 1]
        )
        if min_freq == frequency_indices[list_index - 1]:
            min_freq = freq_idx
    else:
        max_of_prev_band = -1
    # catch if the octave band frequency is equal to another center frequency
    if freq_idx < frequency_indices[-1]:
        if max_freq == frequency_indices[list_index + 1]:
            max_freq = freq_idx
    else:  # if band exceeds cutoff
        max_freq = cutoff_frequency_index
    # catch if the octave band frequency is overlapping with the band below
    while min_freq <= max_of_prev_band:
        min_freq += 1

    band_indices = np.arange(min_freq, max_freq + 1)
    if len(band_indices) == 0:
        band_indices = [freq_idx]

    _logger.debug(f"Band indices from {band_indices[0]} to {band_indices[-1]}")

    return band_indices


def get_positive_frequencies(
    fft_frequencies, fft_coefficients
) -> tuple[np.ndarray, np.ndarray]:
    """Get only positive frequencies and related FFT coefficients.

    Parameters
    ----------
    fft_frequencies: np.ndarray
        Array with frequencies
    fft_coefficients: np.ndarray
        Array with complex-valued FFT coefficients
    """
    signal_length = len(fft_frequencies)
    first_negative_freq = np.argwhere(fft_frequencies < 0)[0, 0]
    frequencies = fft_frequencies[:first_negative_freq]
    fourier_coefficients = fft_coefficients[:first_negative_freq]
    # even signal
    if signal_length % 2 == 0:
        frequencies = np.append(
            frequencies, -1.0 * fft_frequencies[first_negative_freq + 1]
        )
        fourier_coefficients = np.append(
            fourier_coefficients,
            np.conjugate(fft_coefficients[first_negative_freq + 1]),
        )
    return frequencies, fourier_coefficients
