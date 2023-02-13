
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import h5py
import pandas as pd


@dataclass
class ResultImpedance:

    frequency: np.ndarray[float]
    imdedance: np.ndarray[complex]

    def save(self, path: str) -> None:
        data_frame = pd.DataFrame({'frequencies [Hz]': self.frequency,
                                   'impedance [Ohm]': self.imdedance})
        data_frame.to_csv(path, index=False, sep='\t')


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> ResultImpedance:
        pass
