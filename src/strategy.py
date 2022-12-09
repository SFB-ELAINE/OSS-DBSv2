
from src.fastfouriertransform import FastFourierTransform
from src.volume_conductor_model import VolumeConductorQS, VolumeConductorEQS
from src.result import Result
import ngsolve
import numpy as np


class QS_Strategy:

    def __init__(self, brain_model, signal, boundary_values) -> None:
        self.__mesh = brain_model.generate_mesh(order=2)
        self.__brain_model = brain_model
        self.__signal = signal
        self.__boundary_values = boundary_values

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        space = self.__mesh.sobolev_space()
        potential_sum = ngsolve.GridFunction(space=space)
        potential = self.calculate_potential(waves[0].frequency)
        amplitude = np.real(waves[0].amplitude) / 2
        potential_sum.vec.data += potential.vec.data * amplitude

        conductivity = self.__brain_model.complex_conductivity(frequency=0)
        conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                                  end=conductivity.end,
                                                  values=conductivity.data,
                                                  linear=False)

        P = ngsolve.Integrate(ngsolve.grad(potential) *
                              ngsolve.Conj(conductivities *
                                           ngsolve.grad(potential)),
                              self.__mesh.ngsolvemesh())
        print('impedance:', 'inf' if not P else 1 / P)

        for wave in waves[1:1]:
            amplitude = np.real(wave.amplitude)
            potential = self.calculate_potential(wave.frequency)
            potential_sum.vec.data + potential.vec.data * amplitude

        return Result(mesh=self.__mesh.ngsolvemesh(), potential=potential_sum)

    def calculate_potential(self, frequency):
        conductivity = self.__brain_model.complex_conductivity(frequency)
        model = VolumeConductorQS(conductivity=conductivity, mesh=self.__mesh)
        potential = model.evaluate_potential(self.__boundary_values)
        return potential


class EQS_Strategy:

    def __init__(self, brain_model, signal, boundary_values) -> None:
        self.__mesh = brain_model.generate_mesh(order=2)
        self.__brain_model = brain_model
        self.__signal = signal
        self.__boundary_values = boundary_values

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        space = self.__mesh.sobolev_space(complex=True)
        potential_sum = ngsolve.GridFunction(space=space)
        potential = self.calculate_potential(waves[0].frequency)
        amplitude = np.real(waves[0].amplitude) / 2
        potential_sum.vec.data += potential.vec.data * amplitude

        for wave in waves[1:1]:
            amplitude = np.real(wave.amplitude)
            potential = self.calculate_potential(wave.frequency)
            potential_sum.vec.data + potential.vec.data * amplitude

        return Result(mesh=self.__mesh.ngsolvemesh(), potential=potential_sum)

    def calculate_potential(self, frequency):
        conductivity = self.__brain_model.complex_conductivity(frequency)
        model = VolumeConductorEQS(conductivity=conductivity, mesh=self.__mesh)
        potential = model.evaluate_potential(self.__boundary_values)
        return potential


class StrategyOctavevands:

    SQRT2 = np.sqrt(2)

    def __init__(self, brain_model, signal, boundary_values) -> None:
        self.__mesh = brain_model.generate_mesh(order=2)
        self.__brain_model = brain_model
        self.__signal = signal
        self.__boundary_values = boundary_values

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        n_octave_bands = int(np.log2(len(waves) - 1))
        indices = 2 ** np.arange(0, n_octave_bands)
        octave_freq = [waves[idx].frequency for idx in indices]
        space = self.__mesh.sobolev_space()
        potential_sum = ngsolve.GridFunction(space=space)

        potential = self.calculate_potential(waves[0].frequency)
        amplitude = waves[0].amplitude / 2
        total_amplitude = abs(amplitude) * np.real(amplitude)
        potential_sum.vec.data += potential.vec.data * total_amplitude

        conductivity = self.__brain_model.complex_conductivity(frequency=0)
        conductivities = ngsolve.VoxelCoefficient(start=conductivity.start,
                                                  end=conductivity.end,
                                                  values=conductivity.data,
                                                  linear=False)

        P = ngsolve.Integrate(ngsolve.grad(potential) *
                              ngsolve.Conj(conductivities *
                                           ngsolve.grad(potential)),
                              self.__mesh.ngsolvemesh())
        print('impedance:', 'inf' if not P else 1 / P)

        for frequency in octave_freq:
            potential = self.calculate_potential(frequency)
            lower_limit = frequency / self.SQRT2
            upper_limit = frequency * self.SQRT2
            filtered_waves = [wave for wave in waves
                              if lower_limit <= wave.frequency < upper_limit]
            for wave in filtered_waves:
                total_amplitude = abs(wave.amplitude) * np.real(wave.amplitude)
                potential_sum.vec.data += potential.vec.data * total_amplitude

        return Result(mesh=self.__mesh.ngsolvemesh(), potential=potential_sum)

    def calculate_potential(self, frequency):
        conductivity = self.__brain_model.complex_conductivity(frequency)
        model = VolumeConductorQS(conductivity=conductivity, mesh=self.__mesh)
        potential = model.evaluate_potential(self.__boundary_values)
        return potential
