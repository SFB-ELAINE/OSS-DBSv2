
from attr import dataclass
from abc import ABC, abstractmethod
import numpy as np
import ngsolve
import os


@dataclass
class FrequencyComponent:
    fourier_coefficient: float
    frequency: float


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e5

    @abstractmethod
    def result(self, signal, boundary_values, volume_conductor):
        pass

    def _frequency_components(self, signal):
        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        amplitudes = np.fft.rfft(samples) / len(samples) * 2
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        return [FrequencyComponent(fourier_coefficient=amp, frequency=freq)
                for amp, freq in zip(amplitudes, frequencies)]


class NoTruncationTest(SpectrumMode):

    def result(self, signal, boundary_values, volume_conductor):
        freq_components = self._frequency_components(signal)
        frequency = freq_components[77].frequency

        potential, density = volume_conductor.potential(boundary_values,
                                                        frequency)

        amplitude = np.real(freq_components[0].fourier_coefficient) / 2
        potential_sum = potential
        # potential_sum.vec.data = potential.vec.data * amplitude

        for wave in freq_components[1:1]:
            amplitude = np.real(wave.fourier_coefficient)
            potential, density = volume_conductor.evaluate_potential(
                                                        boundary_values,
                                                        wave.frequency)
            potential_sum.vec.data += potential.vec.data * amplitude

        mesh = volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential, density=density)


class Octavevands(SpectrumMode):

    SQRT2 = np.sqrt(2)

    def result(self, signal, boundary_values, volume_conductor):
        freq_components = self._frequency_components(signal)
        indices = 2 ** np.arange(0, int(np.log2(len(freq_components) - 1)))
        octave_freq = [freq_components[idx].frequency for idx in indices]

        potential, density = volume_conductor.potential(boundary_values,
                                                        freq_components[0].frequency)

        amplitude = freq_components[0].fourier_coefficient / 2
        total_amplitude = abs(amplitude) * np.real(amplitude)
        potential_sum = potential
        potential_sum.vec.data += potential.vec.data * total_amplitude

        for frequency in octave_freq:
            potential, P = volume_conductor.evaluate_potential(boundary_values,
                                                               frequency)
            lower_limit = frequency / self.SQRT2
            upper_limit = frequency * self.SQRT2
            amplitudes = [abs(wave.fourier_coefficient)
                          * np.real(wave.fourier_coefficient)
                          for wave in freq_components
                          if lower_limit <= wave.frequency < upper_limit]
            potential_sum.vec.data += potential.vec.data * sum(amplitudes)

        mesh = volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential_sum, density=density)


class Output:

    def __init__(self, mesh, potential, density) -> None:
        self.__mesh = mesh
        self.__potential = potential
        self.__density = density

    def save_mesh(self, path: str = ''):
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

    def save(self, path: str = ''):

        file_base_name = os.path.basename(path)
        file_dir = os.path.dirname(path)

        if not file_base_name:
            file_base_name = 'result'

        if not file_dir:
            file_dir = 'result'

        if not os.path.exists(file_dir):
            os.makedirs(file_dir)

        filename = os.path.join(file_dir, file_base_name)

        field = ngsolve.grad(self.__potential)
        Power = ngsolve.Integrate(field * ngsolve.Conj(self.__density),
                                  self.__mesh)

        print(1 / Power)
        ngsolve.VTKOutput(ma=self.__mesh,
                          coefs=[self.__potential.real,
                                 self.__potential.imag,
                                 field.real,
                                 field.imag,
                                 self.__density.real,
                                 self.__density.imag],
                          names=["potential_real",
                                 "potential_imag",
                                 "field_real",
                                 "field_imag",
                                 "current_density_real",
                                 "current_density_imag"],
                          filename=filename,
                          subdivision=0
                          ).Do()
