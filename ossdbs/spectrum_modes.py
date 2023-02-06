
from attr import dataclass
from abc import ABC, abstractmethod
import numpy as np
import ngsolve
import os
import pickle
from ossdbs.region import Region

@dataclass
class FrequencyComponent:
    fourier_coefficient: float
    frequency: float


class Output:

    def __init__(self, mesh, potential) -> None:
        self.__mesh = mesh
        self.__potential = potential

    def save_mesh(self, path: str = '') -> None:
        file_base_name = os.path.basename(path)
        file_dir = os.path.dirname(path)

        if not file_base_name:
            file_base_name = 'mesh.vol'

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        filename = os.path.join(file_dir, file_base_name)
        self.__mesh.ngmesh.Save(filename)

    def save(self, path: str = '') -> None:

        file_base_name = os.path.basename(path)
        file_dir = os.path.dirname(path)

        if not file_base_name:
            file_base_name = 'result'

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        filename = os.path.join(file_dir, 'potential.data')
        pickler = pickle.Pickler(open(filename, 'wb'))
        pickler.dump([self.__potential, self.__mesh])

        filename = os.path.join(file_dir, file_base_name)
        ngsolve.VTKOutput(ma=self.__mesh,
                          coefs=[self.__potential.real, self.__potential.imag],
                          names=["potential_real", "potential_imag"],
                          filename=filename,
                          subdivision=0
                          ).Do()


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e5

    @abstractmethod
    def result(self, signal, volume_conductor, points) -> Output:
        pass

    def _frequency_components(self, signal):
        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        amplitudes = np.fft.rfft(samples) / len(samples) * 2
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        return [FrequencyComponent(fourier_coefficient=amp, frequency=freq)
                for amp, freq in zip(amplitudes, frequencies)]


class NoTruncationTest(SpectrumMode):

    SQRT2 = np.sqrt(2)

    def result(self, signal, volume_conductor, points):

        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        octave_indices = 2 ** np.arange(0, int(np.log2(len(frequencies) - 1)))
        octave_frequencies = frequencies[octave_indices]
        lower_limits = octave_frequencies / self.SQRT2
        upper_limits = octave_frequencies * self.SQRT2
        starts = lower_limits // signal.frequency + 1
        ends = upper_limits // signal.frequency + 1

        starts = np.concatenate([[0], starts, ends[-1:]]).astype(int)
        ends = np.concatenate([[1], ends, [len(frequencies)]]).astype(int)
        octave_frequencies = np.concatenate([[0],
                                             octave_frequencies,
                                        [octave_frequencies[-1] * self.SQRT2]])

        mesh = volume_conductor.mesh.ngsolvemesh()
        mesh_integrated_points = [mesh(*point) for point in points]
        data = np.zeros((len(points), len(frequencies)))

        # for frequency, start, end in zip(octave_frequencies, starts, ends):
        #     result = volume_conductor.potential(frequency)
        #     potential = result.gridfunction
        #     for index, mip in enumerate(mesh_integrated_points):
        #         value = potential(mip)
        #         data[index, start:end] = value * complex_values[start:end]

        result = volume_conductor.potential(frequencies[77])
        potential = result.gridfunction

        for index, mip in enumerate(mesh_integrated_points):
            value = potential(mip)
            print(value)

        return Output(mesh=mesh,
                      potential=result.gridfunction)


class Octavevands(SpectrumMode):

    SQRT2 = np.sqrt(2)

    def result(self, signal, volume_conductor, points):
        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        complex_values = np.fft.rfft(samples)
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        octave_indices = 2 ** np.arange(0, int(np.log2(len(frequencies) - 1)))
        octave_frequencies = np.append(0, frequencies[octave_indices])
        lower_limits = octave_frequencies / self.SQRT2
        upper_limits = octave_frequencies * self.SQRT2

        frequency = frequencies[0]
        lower_limit = lower_limits[0]
        upper_limit = upper_limits[0]

        # for frequency in octave_frequencies:
        potential = volume_conductor.potential(frequency=frequency)
        fft_values = [value for value, freq in zip(complex_values, frequencies)
                      if lower_limit <= freq < upper_limit]

        


        freq_components = self._frequency_components(signal)

        result = volume_conductor.potential(boundary_values,
                                            freq_components[0].frequency)
        potential, density, impedance = result
        amplitude = freq_components[0].fourier_coefficient / 2
        total_amplitude = abs(amplitude) * np.real(amplitude)
        potential_sum = potential
        potential_sum.vec.data += potential.vec.data * total_amplitude
        impedances = [impedance]
        frequencies = [freq_components[0].frequency]

        for frequency in octave_frequencies:
            result = volume_conductor.evaluate_potential(boundary_values,
                                                         frequency)
            potential, density, impedance = result
            lower_limit = frequency / self.SQRT2
            upper_limit = frequency * self.SQRT2
            
            complex_values = [abs(wave.fourier_coefficient)
                          * np.real(wave.fourier_coefficient)
                          for wave in freq_components
                          if lower_limit <= wave.frequency < upper_limit]
            
            potential_sum.vec.data += potential.vec.data * sum(complex_values)

        mesh = volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh,
                      potential=potential,
                      density=density,
                      impedances=impedances,
                      frequencies=frequencies)
