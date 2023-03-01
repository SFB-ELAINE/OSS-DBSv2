
from ossdbs.stimulation_signal import Signal
from ossdbs.volume_conductor import VolumeConductor
from abc import ABC, abstractmethod
from dataclasses import dataclass
import ngsolve
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

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> Impedances:
        pass


class LogarithmScanning(SpectrumMode):

    def compute(self, signal: Signal, volume_conductor: VolumeConductor,):
        frequencies = self.__frequencies(signal)
        mesh = volume_conductor.mesh.ngsolvemesh()
        impedance = np.zeros(len(frequencies), dtype=complex)

        for index, frequency in enumerate(frequencies[:2]):
            solution = volume_conductor.compute_solution(frequency)
            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            impedance[index] = 1 / power if power else 0
            print(impedance[index])

        return Impedances(frequency=frequencies, imdedance=impedance)

    def __frequencies(self, signal: Signal):
        n_fft_frequencies = len(signal.fft_analysis()[1])
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

        def lower_limit(self):
            return self.frequency / self.SQRT2

        def upper_limit(self):
            return self.frequency * self.SQRT2

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor) -> Impedances:

        frequencies = signal.fft_analysis()[1]
        octave_bands = self.__octave_bands(signal, len(frequencies))
        impedances = np.zeros((len(frequencies)))
        impedances[0] = self.__compute_impedance(volume_conductor, 0)

        for octave_band in octave_bands:
            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            impedance = self.__compute_impedance(volume_conductor,
                                                 octave_band.frequency)
            impedances[start:end] = impedance

        return Impedances(frequency=frequencies, imdedance=impedances)

    def __octave_bands(self, signal, n_frequencies):
        n_octaves = int(np.log2(n_frequencies - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = signal.frequency * octave_indices
        return [self.OctaveBand(freq) for freq in octave_frequencies]

    @staticmethod
    def __compute_impedance(volume_conductor, frequency):
        solution = volume_conductor.compute_solution(frequency)
        field = ngsolve.grad(solution.potential)
        curr_dens_conj = ngsolve.Conj(solution.current_density)
        mesh = volume_conductor.mesh.ngsolvemesh()
        power = ngsolve.Integrate(field * curr_dens_conj, mesh)
        voltage = 1
        return voltage / power if power else 0
