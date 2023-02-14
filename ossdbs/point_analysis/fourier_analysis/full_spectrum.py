
from ossdbs.point_analysis.fourier_analysis.spectrum import SpectrumMode
from ossdbs.point_analysis.fourier_analysis.spectrum import TimeResult
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.stimulation_signal import Signal
import numpy as np


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
