

from ossdbs.electrodes.contacts import Contacts
from ossdbs.point_analysis.field_solution import FieldSolution
from ossdbs.point_analysis.spectrum_modes.spectrum_mode import SpectrumMode
from ossdbs.point_analysis.time_results import TimeResult
from ossdbs.stimmulation_signals.trapzoid_signal import Signal
from ossdbs.fem import VolumeConductor

import numpy as np
import ngsolve


class FullSpectrumBipolar(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                points: np.ndarray,
                contacts: Contacts,
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
            return TimeResult(points=points,
                              potential=np.zeros(len(points)),
                              current_density=np.zeros((len(points), 3)),
                              time_steps=np.array([0]),
                              field_solution=None)

        for index, frequency in enumerate(frequencies):
            solution = volume_conductor.compute_solution(frequency, contacts)
            potential_mip = [solution.potential(mip) for mip in mips]
            current_dens_mip = [solution.current_density(mip) for mip in mips]
            conductivity_mip = [solution.conductivity(mip) for mip in mips]

            if frequency == signal.frequency:
                field_solution = FieldSolution(solution=solution, mesh=ng_mesh)

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
                          time_steps=time_steps,
                          field_solution=field_solution)

    @staticmethod
    def __ifft(fft_spectrum: np.ndarray) -> np.ndarray:
        # inverse fft for only 1000 spectrums at a time
        # to reduce memory stress
        step = 1000
        n_points = fft_spectrum.shape[0]
        return np.concatenate([np.fft.irfft(fft_spectrum[idx:idx + step], axis=1)
                               for idx in range(0, n_points, step)])
