

from ossdbs.electrodes.contacts import Contacts
from ossdbs.stimulation_signal import Signal
from ossdbs.spectrum_modes.voltage_setting import VoltageSetter
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import ngsolve
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
            file.create_dataset(
                'Current_density[A/m2]', data=self.current_density)

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
class FFTResult:

    points: np.ndarray
    frequency: np.ndarray
    potential: np.ndarray
    current_density: np.ndarray
    conductivity: np.ndarray

    def save(self, path: str) -> None:
        with h5py.File(path, "w") as file:
            file.create_dataset("points", data=self.points)
            file.create_dataset("frequencies", data=self.frequency)
            file.create_dataset("potential", data=self.potential)
            file.create_dataset("current_density", data=self.current_density)
            file.create_dataset("conductivity", data=self.conductivity)


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e4

    def __init__(self, voltage_setter: VoltageSetter) -> None:
        self._voltage_setting = voltage_setter

    @abstractmethod
    def compute(self, signal, volume_conductor, points) -> None:
        pass


class FullSpectrum(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                points: np.ndarray,
                contacts: Contacts,
                output: str
                ) -> TimeResult:

        complex_values = signal.fft_analysis()
        frequencies = signal.fft_frequncies()

        source_index = np.argmax([contact.current
                                  for contact in contacts.active()])
        current_set = contacts.active()[source_index].current

        ng_mesh = volume_conductor.mesh.ngsolvemesh()
        included_index = volume_conductor.mesh.is_included(points)
        mips = [ng_mesh(*point) for point in points[included_index]]

        data_shape = len(points), len(frequencies)
        potentials_fft = np.zeros(data_shape, dtype=complex)
        current_dens_fft = np.zeros((*data_shape, 3), dtype=complex)
        conductivities = np.zeros(data_shape, dtype=complex)

        if not len(mips):
            potentials_t = self.__ifft(potentials_fft)
            current_densitys_t = self.__ifft(current_dens_fft)
            sample_spacing = 1 / (signal.frequency * len(potentials_t))
            time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

            return TimeResult(points=points,
                              potential=potentials_t,
                              current_density=current_densitys_t,
                              time_steps=time_steps)

        for index, frequency in enumerate(frequencies[:2]):
            solution = volume_conductor.compute_solution(frequency, contacts)
            potential_mip = [solution.potential(mip) for mip in mips]
            current_dens_mip = [solution.current_density(mip) for mip in mips]
            conductivity_mip = [solution.conductivity(mip) for mip in mips]

            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            mesh = volume_conductor.mesh.ngsolvemesh()
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            voltage = 1
            impedance = voltage / power if power else 0

            pointer = complex_values[index] * impedance * current_set
            pt_index = (included_index, index)
            potentials_fft[pt_index] = np.array(potential_mip) * pointer
            current_dens_fft[pt_index] = np.array(current_dens_mip) * pointer
            conductivities[pt_index] = conductivity_mip

        potentials_t = self.__ifft(potentials_fft)
        current_densitys_t = self.__ifft(current_dens_fft)

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

        return TimeResult(points=points,
                          potential=potentials_t,
                          current_density=current_densitys_t,
                          time_steps=time_steps)

    @staticmethod
    def __ifft(fft_spectrum: np.ndarray) -> np.ndarray:
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        step = 1000
        n_points = fft_spectrum.shape[0]
        return np.concatenate([np.fft.irfft(fft_spectrum[idx:idx+step], axis=1)
                               for idx in range(0, n_points, step)])


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
                points: np.ndarray,
                contacts: Contacts,
                output: str,
                ) -> TimeResult:

        source_index = np.argmax([contact.current
                                  for contact in contacts.active()])
        current_set = contacts.active()[source_index].current

        complex_values = signal.fft_analysis()
        frequencies = signal.fft_frequncies()
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

        if not len(mips):
            potentials_t = self.__ifft(potentials_fft)
            current_densitys_t = self.__ifft(current_dens_fft)
            sample_spacing = 1 / (signal.frequency * len(potentials_t))
            time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

            return TimeResult(points=points,
                              potential=potentials_t,
                              current_density=current_densitys_t,
                              time_steps=time_steps)

        new_contacts = self._voltage_setting.set_voltages(0.0,
                                                          contacts,
                                                          volume_conductor)
        solution = volume_conductor.compute_solution(0.0, new_contacts)
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
            new_contacts = self._voltage_setting.set_voltages(frequency,
                                                              contacts,
                                                              volume_conductor)
            solution = volume_conductor.compute_solution(frequency,
                                                         new_contacts)
            potential_mip = [solution.potential(mip) for mip in mips]
            current_dens_mip = [solution.current_density(mip) for mip in mips]
            conductivity = [solution.conductivity(mip) for mip in mips]

            field = ngsolve.grad(solution.potential)
            curr_dens_conj = ngsolve.Conj(solution.current_density)
            mesh = volume_conductor.mesh.ngsolvemesh()
            power = ngsolve.Integrate(field * curr_dens_conj, mesh)
            voltage = 1
            impedance = voltage / power if power else 0

            start = int(octave_band.lower_limit() / signal.frequency + 1)
            end = int(octave_band.upper_limit() / signal.frequency + 1)
            end = min(end, len(frequencies))

            for index in range(start, end):
                pointer = complex_values[index] * impedance * current_set
                pt_index = (included_index, index)
                potentials_fft[pt_index] = np.array(potential_mip) * pointer
                current_dens_fft[pt_index] = (np.array(current_dens_mip) *
                                              pointer)
                conductivities[pt_index] = conductivity

        potentials_t = self.__ifft(potentials_fft)
        current_densitys_t = self.__ifft(current_dens_fft)
        sample_spacing = 1 / (signal.frequency * len(potentials_t))
        time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

        return TimeResult(points=points,
                          potential=potentials_t,
                          current_density=current_densitys_t,
                          time_steps=time_steps)

    @staticmethod
    def __ifft(fft_spectrum: np.ndarray) -> np.ndarray:
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        step = 1000
        n_points = fft_spectrum.shape[0]
        return np.concatenate([np.fft.irfft(fft_spectrum[idx:idx+step], axis=1)
                               for idx in range(0, n_points, step)])
