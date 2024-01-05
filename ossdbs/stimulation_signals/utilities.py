import multiprocessing as mp
from copy import deepcopy
from functools import partial
from multiprocessing import sharedctypes
from typing import List, Tuple, Union

import numpy as np
from scipy.fft import ifft


def generate_signal(
    coefficients, harmonics, frequency, dt, shift, timesteps, SigmaApprox=False
):
    """Compute signal in parallel.

    TODO document
    """
    global shared_array
    signal = np.ctypeslib.as_ctypes(np.zeros(timesteps, dtype=float))
    shared_array = sharedctypes.RawArray(signal._type_, signal)

    p = mp.Pool()
    time_ind = np.arange(timesteps)
    _ = p.map(
        partial(
            _gen_signal, coefficients, harmonics, frequency, dt, shift, SigmaApprox
        ),
        time_ind,
    )
    signal = np.ctypeslib.as_array(shared_array)
    p.terminate()

    signal = deepcopy(signal)
    return signal


def _gen_signal(coefficient, harmonics, frequency, dt, shift, SigmaApprox, n):
    tmp = np.ctypeslib.as_array(shared_array)
    if SigmaApprox:
        sigma = np.sinc(harmonics / (harmonics[-1] + 1))
    else:
        sigma = np.ones(harmonics.shape)
    signal = 2.0 * np.sum(
        sigma
        * coefficient
        * np.exp(harmonics * 1j * 2.0 * np.pi * frequency * (n * dt - shift))
    )
    signal -= coefficient[0]
    tmp[n] = np.real(signal)


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
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute time-domain signal via fft."""
    # double the cutoff_frequency to actually sample until there
    signal = ifft(fft_signal).real
    timesteps = get_timesteps(cutoff_frequency, base_frequency, len(signal))
    return timesteps, signal


def get_octave_band_indices(frequencies: np.ndarray) -> np.ndarray:
    """Return indices of octave band frequencies."""
    n_octaves = int(np.log2(len(frequencies) - 1)) + 1
    octave_indices = 2 ** np.arange(0, n_octaves)
    return octave_indices


def get_minimum_octave_band_index(freq_idx: int):
    """Get index of lowest frequency in octave band."""
    return int(np.round(freq_idx / np.sqrt(2)))


def get_maximum_octave_band_index(freq_idx: int):
    """Get index of highest frequency in octave band."""
    return int(np.round(freq_idx * np.sqrt(2)))


def get_indices_in_octave_band(
    freq_idx: int, frequency_indices: list, cutoff_frequency_index: int
) -> Union[List, np.ndarray]:
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
    return band_indices
