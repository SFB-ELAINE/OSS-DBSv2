
from attr import dataclass
from abc import ABC, abstractmethod
import numpy as np
import ngsolve
import os


@dataclass
class FFTWave:
    amplitude: float
    frequency: float


class SpectrumMode(ABC):

    SPACING_FACTOR = 1e5

    @abstractmethod
    def result(self, signal, boundary_values, volume_conductor):
        pass

    def _fft_waves(self, signal):
        sample_spacing = 1 / (signal.frequency * self.SPACING_FACTOR)
        samples = signal.generate_samples(sample_spacing)
        amplitudes = np.fft.rfft(samples) / len(samples) * 2
        frequencies = np.fft.rfftfreq(len(samples), sample_spacing)
        return [FFTWave(amplitude=amp, frequency=freq)
                for amp, freq in zip(amplitudes, frequencies)]


class NoTruncationTest(SpectrumMode):

    def result(self, signal, boundary_values, volume_conductor):
        waves = self._fft_waves(signal)
        frequency = waves[77].frequency

        potential, density = volume_conductor.evaluate_potential(
                                                               boundary_values,
                                                               frequency
                                                               )

        amplitude = np.real(waves[0].amplitude) / 2
        potential_sum = potential
        # potential_sum.vec.data = potential.vec.data * amplitude

        for wave in waves[1:1]:
            amplitude = np.real(wave.amplitude)
            potential, density = volume_conductor.evaluate_potential(
                                                        boundary_values,
                                                        wave.frequency)
            potential_sum.vec.data += potential.vec.data * amplitude

        mesh = volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential, density=density)


class Octavevands(SpectrumMode):

    SQRT2 = np.sqrt(2)

    def result(self, signal, boundary_values, volume_conductor):
        waves = self._fft_waves(signal)
        indices = 2 ** np.arange(0, int(np.log2(len(waves) - 1)))
        octave_freq = [waves[idx].frequency for idx in indices]

        potential, density = volume_conductor.evaluate_potential(
                                                            boundary_values,
                                                            waves[0].frequency)

        amplitude = waves[0].amplitude / 2
        total_amplitude = abs(amplitude) * np.real(amplitude)
        potential_sum = potential
        potential_sum.vec.data += potential.vec.data * total_amplitude

        for frequency in octave_freq:
            potential, P = volume_conductor.evaluate_potential(boundary_values,
                                                               frequency)
            lower_limit = frequency / self.SQRT2
            upper_limit = frequency * self.SQRT2
            amplitudes = [abs(wave.amplitude) * np.real(wave.amplitude)
                          for wave in waves
                          if lower_limit <= wave.frequency < upper_limit]
            potential_sum.vec.data += potential.vec.data * sum(amplitudes)

        mesh = volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential_sum, density=density)


class Output:

    def __init__(self, mesh, potential, density) -> None:
        self.__mesh = mesh
        self.__potential = potential
        self.__density = density

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
