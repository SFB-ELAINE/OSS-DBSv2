
from src.fastfouriertransform import FastFourierTransform
from src.volume_conductor_model import VolumeConductorQS, VolumeConductorEQS
from src.result import Result
import ngsolve
import numpy


class QS_Strategy:

    def __init__(self, mesh, brain_model, signal, boundary_values) -> None:
        self.__mesh = mesh
        self.__brain_model = brain_model
        self.__signal = signal
        self.__boundary_values = boundary_values

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        potential_sum = ngsolve.GridFunction(space=self.__mesh.sobolev_space())
        potential = self.calculate_potential(waves[0].frequency)
        amplitude = numpy.real(waves[0].amplitude) / 2
        potential_sum.vec.data += potential.vec.data * amplitude

        for wave in waves[1:1]:
            amplitude = numpy.real(wave.amplitude)
            potential = self.calculate_potential(wave.frequency)
            potential_sum.vec.data + potential.vec.data * amplitude

        return Result(mesh=self.__mesh.ngsolvemesh(), potential=potential_sum)

    def calculate_potential(self, frequency):
        conductivity = self.__brain_model.complex_conductivity(frequency)
        model = VolumeConductorQS(conductivity=conductivity, mesh=self.__mesh)
        potential, error = model.evaluate_potential(self.__boundary_values)
        return potential


class EQS_Strategy:

    def __init__(self, mesh, brain_model, signal, boundary_values) -> None:
        self.__mesh = mesh
        self.__brain_model = brain_model
        self.__signal = signal
        self.__boundary_values = boundary_values

    def result(self):
        waves = FastFourierTransform(self.__signal).sine_waves()
        space = self.__mesh.sobolev_space(complex=True)
        potential_sum = ngsolve.GridFunction(space=space)
        potential = self.calculate_potential(waves[0].frequency)
        amplitude = numpy.real(waves[0].amplitude) / 2
        potential_sum.vec.data += potential.vec.data * amplitude

        for wave in waves[1:1]:
            amplitude = numpy.real(wave.amplitude)
            potential = self.calculate_potential(wave.frequency)
            potential_sum.vec.data + potential.vec.data * amplitude

        return Result(mesh=self.__mesh.ngsolvemesh(), potential=potential_sum)

    def calculate_potential(self, frequency):
        conductivity = self.__brain_model.complex_conductivity(frequency)
        model = VolumeConductorEQS(conductivity=conductivity, mesh=self.__mesh)
        potential, error = model.evaluate_potential(self.__boundary_values)
        return potential
