from dataclasses import dataclass
from src.brainsubstance import BrainSubstance
from src.dielectric_model.dielectric_model import DielectricModel as Model
import numpy as np


class Parameters:
    pass


@dataclass
class WhiteMatterParameters(Parameters):
    alpha: float = 0.24
    epsilon: float = 37.0
    epsilon_inf: float = 4.0
    sigma: float = 0.47
    tau: float = 8.04e-12


@dataclass
class GrayMatterParameters(Parameters):
    alpha: float = 0.12
    epsilon: float = 55.5
    epsilon_inf: float = 4.0
    sigma: float = 1.03
    tau: float = 7.76e-12


class DielectricModelCSF(Model):

    def permitivity(self, frequency: float) -> float:
        return 80

    def conductivity(self, frequency: float) -> float:
        return 1.79

    @classmethod
    def create_model(cls, material: BrainSubstance) -> 'DielectricModel':

        if material is not BrainSubstance.CEREBROSPINAL_FLUID:
            material_parameters = {
                            BrainSubstance.WHITE_MATTER: WhiteMatterParameters,
                            BrainSubstance.GRAY_MATTER: GrayMatterParameters}

            return cls(material_parameters[material])

        return cls()


class DielectricModel(Model):
    """Model variant 2 for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: Parameters) -> None:
        self.__epsilon_inf = parameters.epsilon_inf
        self.__sigma = parameters.sigma
        self.__alpha = parameters.alpha
        self.__epsilon = parameters.epsilon
        self.__tau = parameters.tau

    def __complex_permitivity(self, frequency: float) -> complex:
        divisor = 1 + (1j * frequency * self.__tau) ** (1 - self.__alpha)
        eps_dispersion = self.__epsilon_inf +\
            (self.__epsilon - self.__epsilon_inf) / divisor

        if frequency == 0:
            return eps_dispersion

        conductivity_term = self.__sigma / (1j * frequency * self.e0)
        return eps_dispersion + conductivity_term

    def permitivity(self, frequency: float) -> float:
        return np.real(self.__complex_permitivity(frequency))

    def conductivity(self, frequency: float) -> float:
        if frequency == 0:
            return self.__sigma

        permitivity_imaginary = np.imag(self.__complex_permitivity(frequency))
        return -1 * frequency * self.e0 * permitivity_imaginary

    @classmethod
    def create_model(cls, material: BrainSubstance = BrainSubstance(1)) \
            -> 'DielectricModel':

        if material is BrainSubstance.CEREBROSPINAL_FLUID:
            return DielectricModelCSF()

        material_parameters = {
                            BrainSubstance.WHITE_MATTER: WhiteMatterParameters,
                            BrainSubstance.GRAY_MATTER: GrayMatterParameters}

        return cls(material_parameters[material])