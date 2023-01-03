
from ossdbs.fastfouriertransform import FastFourierTransform
import numpy as np
import ngsolve
import os


class AllFrequenciesStrategy:

    def __init__(self, signal, boundary_values, volume_conductor) -> None:
        self.__signal = signal
        self.__boundary_values = boundary_values
        self.__volume_conductor = volume_conductor

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        frequency = waves[77].frequency

        potential, density = self.__volume_conductor.evaluate_potential(
                                                        self.__boundary_values,
                                                        frequency)

        amplitude = np.real(waves[0].amplitude) / 2
        potential_sum = potential
        # potential_sum.vec.data = potential.vec.data * amplitude

        for wave in waves[1:1]:
            amplitude = np.real(wave.amplitude)
            potential, density = self.__volume_conductor.evaluate_potential(
                                                        self.__boundary_values,
                                                        wave.frequency)
            potential_sum.vec.data += potential.vec.data * amplitude

        mesh = self.__volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential, density=density)


class StrategyOctavevands:

    SQRT2 = np.sqrt(2)

    def __init__(self, signal, boundary_values, volume_conductor) -> None:
        self.__signal = signal
        self.__boundary_values = boundary_values
        self.__volume_conductor = volume_conductor

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        indices = 2 ** np.arange(0, int(np.log2(len(waves) - 1)))
        octave_freq = [waves[idx].frequency for idx in indices]

        potential, density = self.__volume_conductor.evaluate_potential(
                                                        self.__boundary_values,
                                                        waves[0].frequency)

        amplitude = waves[0].amplitude / 2
        total_amplitude = abs(amplitude) * np.real(amplitude)
        potential_sum = potential
        potential_sum.vec.data += potential.vec.data * total_amplitude

        for frequency in octave_freq:
            potential, P = self.__volume_conductor.evaluate_potential(
                                                        self.__boundary_values,
                                                        frequency)
            lower_limit = frequency / self.SQRT2
            upper_limit = frequency * self.SQRT2
            filtered_waves = [wave for wave in waves
                              if lower_limit <= wave.frequency < upper_limit]
            for wave in filtered_waves:
                total_amplitude = abs(wave.amplitude) * np.real(wave.amplitude)
                potential_sum.vec.data += potential.vec.data * total_amplitude

        mesh = self.__volume_conductor.mesh.ngsolvemesh()
        return Output(mesh=mesh, potential=potential, density=density)


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
