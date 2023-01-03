
from dataclasses import dataclass
from ossdbs.signals import Signal
import numpy as np


@dataclass
class FFTWave:
    amplitude: float
    frequency: float


class FastFourierTransform:

    SPACING_FACTOR = 1e5

    def __init__(self, signal: Signal) -> None:
        self.__signal = signal

    def sine_waves(self) -> list:
        sample_spacing = 1 / (self.__signal.frequency * self.SPACING_FACTOR)
        samples = self.__signal.generate_samples(sample_spacing)
        amplitudes = np.fft.rfft(samples) / len(samples) * 2
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        return [FFTWave(amplitude=amp, frequency=freq)
                for amp, freq in zip(amplitudes, frequencies)]
