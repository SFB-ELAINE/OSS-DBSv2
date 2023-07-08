import numpy as np
import multiprocessing as mp
from multiprocessing import sharedctypes
from functools import partial
from copy import deepcopy


def generate_signal(coefficients, harmonics, frequency, dt, shift, timesteps, SigmaApprox=False):
    """Compute signal in parallel

    TODO document
    """
    global shared_array
    signal = np.ctypeslib.as_ctypes(np.zeros(timesteps, dtype=float))
    shared_array = sharedctypes.RawArray(signal._type_, signal)

    p = mp.Pool()
    time_ind = np.arange(timesteps)
    _ = p.map(partial(_gen_signal, coefficients, harmonics, frequency, dt, shift, SigmaApprox), time_ind)
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
    signal = 2. * np.sum(sigma * coefficient * np.exp(harmonics * 1j * 2. * np.pi * frequency * (n * dt - shift)))
    signal -= coefficient[0]
    tmp[n] = np.real(signal)


def adjust_cutoff_frequency(cutoff_frequency, frequency):
    """Function to make cutoff frequency multiple of stimulation frequency
    """
    return cutoff_frequency - cutoff_frequency % frequency
