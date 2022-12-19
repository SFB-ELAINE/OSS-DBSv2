
from src.fastfouriertransform import FastFourierTransform
import numpy as np

from src.output import OutputTest


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

        mesh = self.__volume_conductor.mesh
        return OutputTest(mesh=mesh.ngsolvemesh(),
                          potential=potential,
                          density=density)


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

        mesh = self.__volume_conductor.mesh
        return OutputTest(mesh=mesh.ngsolvemesh(),
                          potential=potential,
                          density=density)
