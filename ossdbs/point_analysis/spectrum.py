

from ossdbs.stimulation_signal import Signal
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import h5py


@dataclass
class TimeResult:
    points: np.ndarray
    time_steps: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset('TimeSteps[s]', data=self.time_steps)
            file.create_dataset('Points[mm]', data=self.points)
            file.create_dataset('Potential[V]', data=self.potential)
            file.create_dataset('Current_density[A/m2]',
                                data=self.current_density)

    def save_by_categories(self, path: str, categories: list) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset('TimeSteps', data=self.time_steps)
            start = 0
            for category in categories:
                name, n_points = category
                end = start + n_points
                h5_group = file.create_group(name)
                points = self.points[start:end]
                h5_group.create_dataset('Points', data=points)
                potential = self.potential[start:end]
                h5_group.create_dataset('Potential', data=potential)
                current_density = self.current_density[start:end]
                h5_group.create_dataset('CurrentDensity', data=current_density)


@dataclass
class Result:

    points: np.ndarray[float]
    frequency: np.ndarray[float]
    potential: np.ndarray[complex]
    current_density: np.ndarray[complex]
    conductivity: np.ndarray[complex]

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset("points", data=self.points)
            file.create_dataset("frequencies", data=self.frequency)
            file.create_dataset("potential", data=self.potential)
            file.create_dataset("current_density", data=self.current_density)
            file.create_dataset("conductivity", data=self.conductivity)


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> None:
        pass



class FullSpectrum(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                points: np.ndarray):

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)

        ng_mesh = volume_conductor.mesh.ngsolvemesh()
        included_index = volume_conductor.mesh.is_included(points)
        mips = [ng_mesh(*point) for point in points[included_index]]

        data_shape = len(points), len(frequencies)
        potentials_fft = np.zeros(data_shape, dtype=complex)
        current_dens_fft = np.zeros((*data_shape, 3), dtype=complex)
        conductivities = np.zeros(data_shape, dtype=complex)

        for index, frequency in enumerate(frequencies[:0]):
            solution = volume_conductor.compute_solution(frequency)
            potential_mip = [solution.potential(mip) for mip in mips]
            current_dens_mip = [solution.current_density(mip) for mip in mips]
            conductivity_mip = [solution.conductivity(mip) for mip in mips]

            pointer = complex_values[index]
            pt_index = (included_index, index)
            potentials_fft[pt_index] = np.array(potential_mip) * pointer
            current_dens_fft[pt_index] = np.array(current_dens_mip) * pointer
            conductivities[pt_index] = conductivity_mip

        potentials_t = self.__ifft(potentials_fft)
        current_densitys_t = self.__ifft(current_dens_fft)
        time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

        return TimeResult(points=points,
                          potential=potentials_t,
                          current_density=current_densitys_t,
                          time_steps=time_steps)

    def __ifft(self, fft_spectrum: np.ndarray) -> np.ndarray:
        time_signals = []
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        for start in range(0, fft_spectrum.shape[0], 1000):
            time_signal = np.fft.irfft(fft_spectrum[start:start+1000], axis=1)
            time_signals.append(time_signal)

        return np.concatenate(time_signals, axis=0)


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
                volume_conductor: VolumeConductor,
                points: np.ndarray):

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = signal.frequency * octave_indices
        octave_bands = [self.OctaveBand(freq) for freq in octave_frequencies]

        ng_mesh = volume_conductor.mesh.ngsolvemesh()
        included_index = volume_conductor.mesh.is_included(points)
        mips = [ng_mesh(*point) for point in points[included_index]]

        data_shape = len(points), len(frequencies)
        potentials_fft = np.zeros(data_shape, dtype=complex)
        current_dens_fft = np.zeros((*data_shape, 3), dtype=complex)
        conductivities = np.zeros(data_shape, dtype=complex)

        solution = volume_conductor.compute_solution(0)
        potential_mip = [solution.potential(mip) for mip in mips]
        current_dens_mip = [solution.current_density(mip) for mip in mips]
        conductivity = [solution.conductivity(mip) for mip in mips]

        pt_index = (included_index, 0)
        pointer = complex_values[0]
        potentials_fft[pt_index] = np.array(potential_mip) * pointer
        current_dens_fft[pt_index] = np.array(current_dens_mip) * pointer
        conductivities[pt_index] = conductivity

        for octave_band in octave_bands:
            frequency = octave_band.frequency
            solution = volume_conductor.compute_solution(frequency)
            potential_mip = [solution.potential(mip) for mip in mips]
            current_dens_mip = [solution.current_density(mip) for mip in mips]
            conductivity = [solution.conductivity(mip) for mip in mips]

            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            end = min(end, len(frequencies))

            for index in range(start, end):
                pointer = complex_values[index]
                pt_index = (included_index, index)
                potentials_fft[pt_index] = np.array(potential_mip) * pointer
                current_dens_fft[pt_index] = (np.array(current_dens_mip) *
                                              pointer)
                conductivities[pt_index] = conductivity

        potentials_t = self.__ifft(potentials_fft)
        current_densitys_t = self.__ifft(current_dens_fft)
        time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

        return TimeResult(points=points,
                          potential=potentials_t,
                          current_density=current_densitys_t,
                          time_steps=time_steps)

    def __ifft(self, fft_spectrum: np.ndarray) -> np.ndarray:
        time_signals = []
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        for start in range(0, fft_spectrum.shape[0], 1000):
            time_signal = np.fft.irfft(fft_spectrum[start:start+1000], axis=1)
            time_signals.append(time_signal)

        return np.concatenate(time_signals, axis=0)
