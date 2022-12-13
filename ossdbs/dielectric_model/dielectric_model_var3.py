from dataclasses import dataclass
from ossdbs.brainsubstance import Material
from ossdbs.dielectric_model.dielectric_model import AbstractDielectricModel
import numpy as np


class Paramters:
    pass


@dataclass
class WhiteMatterParameters(Paramters):
    alpha_1: float = 0.46
    alpha_2: float = 0.05
    eps_delta_1: float = 4.29e3
    eps_delta_2: float = 32.62e6
    eps_inf: float = 10.27
    sigma: float = 0.0192
    tau_1: float = 3.71e-6
    tau_2: float = 8.29e-3


@dataclass
class GrayMatterParameters(Paramters):
    alpha_1: float = 0.45
    alpha_2: float = 0.06
    eps_delta_1: float = 3.68e3
    eps_delta_2: float = 42.26e6
    eps_inf: float = 1.0
    sigma: float = 0.027
    tau_1: float = 1.14e-6
    tau_2: float = 5.83e-3


class DielectricModelCSF(AbstractDielectricModel):

    def relative_permitivity(self, frequency: float) -> float:
        return 80

    def conductivity(self, frequency: float) -> float:
        return 1.79

    @classmethod
    def create_model(cls, material: Material) -> 'DielectricModel':

        if material is not Material.CSF:
            material_parameters = {
                            Material.WHITE_MATTER: WhiteMatterParameters,
                            Material.GRAY_MATTER: GrayMatterParameters}

            return cls(material_parameters[material])

        return cls()


class DielectricModel(AbstractDielectricModel):
    """"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: Paramters) -> None:
        self.__eps_inf = parameters.eps_inf
        self.__sigma = parameters.sigma
        self.__alpha = np.array([parameters.alpha_1,
                                 parameters.alpha_2])
        self.__eps_delta = np.array([parameters.eps_delta_1,
                                     parameters.eps_delta_2])
        self.__tau = np.array([parameters.tau_1,
                               parameters.tau_2])

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
    def create_model(cls, material: Material = Material(1)) \
            -> 'DielectricModel':

        if material is Material.CSF:
            return DielectricModelCSF()

        material_parameters = {
                            Material.WHITE_MATTER: WhiteMatterParameters,
                            Material.GRAY_MATTER: GrayMatterParameters}

        return cls(material_parameters[material])
