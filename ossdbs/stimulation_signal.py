
from ossdbs.signals import Signal
import numpy as np


class StimulationSiganl:

    def __init__(self,
                 signal: Signal,
                 counter_signal: Signal,
                 inter_pulse_width: float) -> None:
        self.__signal = signal
        self.__counter_signal = counter_signal
        self.__inter_pulse_width = inter_pulse_width

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        samples = self.__signal.generate_samples(sample_spacing=sample_spacing)
        