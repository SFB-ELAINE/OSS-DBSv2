
from ossdbs.electrodes.contacts import Contacts
from ossdbs.point_analysis.field_solution import FieldSolution
from ossdbs.point_analysis.spectrum_modes.spectrum_mode import SpectrumMode
from ossdbs.stimmulation_signals.trapzoid_signal import Signal
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.fem import VolumeConductor

import numpy as np


class OctaveBand(SpectrumMode):

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
                ) -> TimeResult:

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
            return TimeResult(points=points,
                              potential=np.zeros(len(points)),
                              current_density=np.zeros((len(points), 3)),
                              time_steps=np.array([0]),
                              field_solution=None)

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

            if frequency == signal.frequency:
                field_solution = FieldSolution(solution=solution, mesh=ng_mesh)

        potentials_t = self.__ifft(potentials_fft)
        current_densitys_t = self.__ifft(current_dens_fft)
        sample_spacing = 1 / (signal.frequency * len(potentials_t))
        time_steps = np.arange(self.SPACING_FACTOR) * sample_spacing

        return TimeResult(points=points,
                          potential=potentials_t,
                          current_density=current_densitys_t,
                          time_steps=time_steps,
                          field_solution=field_solution)

    @staticmethod
    def __ifft(fft_spectrum: np.ndarray) -> np.ndarray:
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        step = 1000
        n_points = fft_spectrum.shape[0]
        return np.concatenate([np.fft.irfft(fft_spectrum[idx:idx+step], axis=1)
                               for idx in range(0, n_points, step)])
