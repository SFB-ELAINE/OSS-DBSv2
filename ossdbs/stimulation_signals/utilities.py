import multiprocessing as mp
from copy import deepcopy
from functools import partial
from multiprocessing import sharedctypes
from scipy.fft import ifft

import numpy as np


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


def retrieve_time_domain_signal_from_fft(
    fft_signal: np.ndarray, cutoff_frequency: float, base_frequency: float
) -> tuple[np.ndarray, np.ndarray]:
    """Compute time-domain signal via fft
    """
    # double the cutoff_frequency to actually sample until there
    cutoff_frequency = adjust_cutoff_frequency(2.0 * cutoff_frequency, base_frequency)
    dt = 1.0 / cutoff_frequency
    signal = ifft(fft_signal).real
    timesteps = dt * np.arange(len(signal))
    return timesteps, signal
