from dataclasses import dataclass
from src.brainsubstance import BrainSubstance
from src.dielectric_model.dielectric_model import DielectricModel as Model
import numpy as np


class Paramters:
    pass


@dataclass
class WhiteMatterParameters(Paramters):
    alpha_1: float = 0.1
    alpha_2: float = 0.1
    alpha_3: float = 0.3
    alpha_4: float = 0.02
    eps_delta_1: float = 32.0
    eps_delta_2: float = 100.0
    eps_delta_3: float = 4.0e4
    eps_delta_4: float = 3.5e7
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau_1: float = 7.96e-12
    tau_2: float = 7.96e-9
    tau_3: float = 53.05e-6
    tau_4: float = 7.958e-3


@dataclass
class GrayMatterParameters(Paramters):
    alpha_1: float = 0.1
    alpha_2: float = 0.1
    alpha_3: float = 0.22
    alpha_4: float = 0.0
    eps_delta_1: float = 45.0
    eps_delta_2: float = 400.0
    eps_delta_3: float = 2.0e5
    eps_delta_4: float = 4.5e7
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau_1: float = 7.96e-12
    tau_2: float = 15.92e-9
    tau_3: float = 106.1e-6
    tau_4: float = 5.305e-3


class DielectricModelCSF(Model):

    def relative_permitivity(self, frequency: float) -> float:
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
    """Model variant 1 for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: Paramters) -> None:
        self.__eps_inf = parameters.eps_inf
        self.__sigma = parameters.sigma
        self.__alpha = np.array([parameters.alpha_1,
                                 parameters.alpha_2,
                                 parameters.alpha_3,
                                 parameters.alpha_4])
        self.__eps_delta = np.array([parameters.eps_delta_1,
                                     parameters.eps_delta_2,
                                     parameters.eps_delta_3,
                                     parameters.eps_delta_4])
        self.__tau = np.array([parameters.tau_1,
                               parameters.tau_2,
                               parameters.tau_3,
                               parameters.tau_4])

    def __complex_permitivity(self, frequency: float) -> complex:
        divisor = 1 + (1j * frequency * self.__tau) ** (1 - self.__alpha)
        eps_dispersion = self.__eps_inf + np.sum(self.__eps_delta / divisor)

        if frequency == 0:
            return eps_dispersion

        conductivity_term = self.__sigma / (1j * frequency * self.e0)
        return eps_dispersion + conductivity_term

    def relative_permitivity(self, frequency: float) -> float:
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
