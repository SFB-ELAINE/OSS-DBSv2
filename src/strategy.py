
from src.fastfouriertransform import FastFourierTransform
import numpy as np


class AllFrequenciesStrategy:

    def __init__(self, signal, boundary_values, volume_conductor) -> None:
        self.__signal = signal
        self.__boundary_values = boundary_values
        self.__volume_conductor = volume_conductor

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        potential, P = self.__volume_conductor.evaluate_potential(
            self.__boundary_values,
            waves[0].frequency)
        amplitude = np.real(waves[0].amplitude) / 2
        potential_sum = potential
        potential_sum.vec.data = potential.vec.data * amplitude

        print('impedance:', 'inf' if not P else 1 / P)

        for wave in waves[1:1]:
            amplitude = np.real(wave.amplitude)
            potential, _ = self.__volume_conductor.evaluate_potential(
                self.__boundary_values,
                wave.frequency)
            potential_sum.vec.data + potential.vec.data * amplitude

        return potential_sum


class StrategyOctavevands:

    SQRT2 = np.sqrt(2)

    def __init__(self, signal, boundary_values, volume_conductor) -> None:
        self.__signal = signal
        self.__boundary_values = boundary_values
        self.__volume_conductor = volume_conductor

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        n_octave_bands = int(np.log2(len(waves) - 1))
        indices = 2 ** np.arange(0, n_octave_bands)
        octave_freq = [waves[idx].frequency for idx in indices]

        potential, P = self.__volume_conductor.evaluate_potential(
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

        return potential_sum
