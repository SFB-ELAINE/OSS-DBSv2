
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class Impedances:

    frequency: np.ndarray
    imdedance: np.ndarray

    def save(self, path: str) -> None:
        data_frame = pd.DataFrame({'frequencies [Hz]': self.frequency,
                                   'Resistance [Ohm]': np.real(self.imdedance),
                                   'Reactance [Ohm]': np.imag(self.imdedance)})
        data_frame.to_csv(path, index=False, sep=',')


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> Impedances:
        pass
