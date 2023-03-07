
from typing import List
from ossdbs.contacts import Contacts
from ossdbs.stimulation_signal import Signal
from ossdbs.volume_conductor import VolumeConductor
from abc import ABC, abstractmethod
from dataclasses import dataclass
import ngsolve
import numpy as np
import pandas as pd


@dataclass
class Impedances:

    frequencies: np.ndarray
    imdedances: np.ndarray

    def save(self, path: str) -> None:
        data = pd.DataFrame({'frequencies [Hz]': self.frequencies,
                             'Resistance [Ohm]': np.real(self.imdedances),
                             'Reactance [Ohm]': np.imag(self.imdedances)}
                            )
        data.to_csv(path, index=False, sep=',')


class SpectrumMode(ABC):

    @abstractmethod
    def compute(self, signal, volume_conductor, contacts) -> Impedances:
        pass


class LogarithmScanning(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                contacts: Contacts
                ) -> Impedances:
        frequencies = self.__frequencies(signal)
        mesh = volume_conductor.mesh.ngsolvemesh()
        impedance = np.zeros(len(frequencies), dtype=complex)

        for index, frequency in enumerate(frequencies[:2]):
            solution = volume_conductor.compute_solution(frequency, contacts)
            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            impedance[index] = 1 / power if power else 0
            print(impedance[index])

        return Impedances(frequencies=frequencies, imdedances=impedance)

    def __frequencies(self, signal: Signal) -> List[float]:
        n_fft_frequencies = len(signal.fft_frequncies())
        n_intervalls = int(np.log10(n_fft_frequencies)) + 1
        frequencies = [0.0]
        for index in range(n_intervalls):
            interval__freq = signal.frequency * 10 ** index * np.arange(1, 10)
            frequencies.extend(interval__freq)
        return frequencies


class OctaveBandMode(SpectrumMode):

    class OctaveBand:

        SQRT2 = np.sqrt(2)

        def __init__(self, frequency: float) -> None:
            self.frequency = frequency

        def lower_limit(self) -> float:
            return self.frequency / self.SQRT2

        def upper_limit(self) -> float:
            return self.frequency * self.SQRT2

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                contacts: Contacts
                ) -> Impedances:

        frequencies = signal.fft_frequncies()
        impedances = np.zeros((len(frequencies)))
        impedances[0] = self.__compute_impedance(volume_conductor, 0, contacts)

        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        for octave_band in self.__octave_bands(signal.frequency, n_octaves):
            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            impedance = self.__compute_impedance(volume_conductor,
                                                 octave_band.frequency)
            impedances[start:end] = impedance

        return Impedances(frequencies=frequencies, imdedances=impedances)

    def __octave_bands(self,
                       frequency: float,
                       n_octaves: int
                       ) -> List[OctaveBand]:
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = frequency * octave_indices
        return [self.OctaveBand(freq) for freq in octave_frequencies]

    @staticmethod
    def __compute_impedance(volume_conductor: VolumeConductor,
                            frequency: float,
                            contacts: Contacts
                            ) -> float:
        solution = volume_conductor.compute_solution(frequency, contacts)
        field = ngsolve.grad(solution.potential)
        curr_dens_conj = ngsolve.Conj(solution.current_density)
        mesh = volume_conductor.mesh.ngsolvemesh()
        power = ngsolve.Integrate(field * curr_dens_conj, mesh)
        voltage = 1
        return voltage / power if power else 0
