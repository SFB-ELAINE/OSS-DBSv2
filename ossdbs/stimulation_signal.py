
from ossdbs.signals import Signal
import numpy as np


class StimulationSiganl:

    def __init__(self,
                 signal: Signal,
                 counter: Signal,
                 counter_offset: float) -> None:
        self.__signal = signal
        self.__counter = counter
        self.__counter_offset = counter_offset

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        samples = self.__signal.generate_samples(sample_spacing=sample_spacing)
        n_samples = len(samples)
        offset_samples = n_samples * self.__counter_offset
        counter_samples = -1 * self.__counter.generate_samples(sample_spacing)
        samples[offset_samples:] = counter_samples[-offset_samples]
        return samples
