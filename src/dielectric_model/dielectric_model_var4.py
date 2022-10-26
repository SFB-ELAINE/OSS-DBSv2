from dataclasses import dataclass
from src.brainsubstance import BrainSubstance
from src.dielectric_model.dielectric_model import DielectricModel
import numpy as np


class Paramters:
    pass


@dataclass
class WhiteMatterParameters(Paramters):
    alpha: float = 0.02
    eps_delta: float = 49.03e3
    eps_inf: float = 7.77
    sigma: float = 0.0611
    tau: float = 159.77e-6


@dataclass
class GrayMatterParameters(Paramters):
    alpha: float = 0.0
    eps_delta: float = 95.55e3
    eps_inf: float = 11.92
    sigma: float = 0.107
    tau: float = 170.91e-6


class DielectricModelVariant1(DielectricModel):
    """Model variant 1 for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: Paramters) -> None:
        self.__eps_inf = parameters.eps_inf
        self.__sigma = parameters.sigma
        self.__alpha = parameters.alpha
        self.__eps_delta = parameters.eps_delta
        self.__tau = parameters.tau

    def __complex_permitivity(self, frequency: float) -> complex:
        divisor = 1 + (1j * frequency * self.__tau) ** (1 - self.__alpha)
        eps_dispersion = self.__eps_inf + np.sum(self.__eps_delta / divisor)

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
        material_parameters = {
                            BrainSubstance.WHITE_MATTER: WhiteMatterParameters,
                            BrainSubstance.GRAY_MATTER: GrayMatterParameters}

        return cls(material_parameters[material])
