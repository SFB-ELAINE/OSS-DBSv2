from dataclasses import dataclass
from ossdbs.materials import Material
from ossdbs.dielectric_model import DielectricModel
from abc import ABC
import numpy as np


class ColeColeParamters(ABC):
    alpha: np.ndarray
    eps_delta: np.ndarray
    eps_inf: float
    sigma: float
    tau: np.ndarray


@dataclass
class WhiteMatterParameters(ColeColeParamters):
    alpha: np.ndarray = np.array([0.1, 0.1, 0.3, 0.02])
    eps_delta: np.ndarray = np.array([32.0, 100.0, 4.0e4, 3.5e7])
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau: np.ndarray = np.array([7.958e-12, 7.958e-9, 53.052e-6, 7.958e-3])


@dataclass
class GrayMatterParameters(ColeColeParamters):
    alpha: np.ndarray = np.array([0.1, 0.15, 0.22, 0.0])
    eps_delta: np.ndarray = np.array([45.0, 400.0, 2.0e5, 4.5e7])
    eps_inf: float = 4.0
    sigma: float = 0.02
    tau: np.ndarray = np.array([7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3])


@dataclass
class CerebroSpinalFluidParameters(ColeColeParamters):
    alpha: np.ndarray = np.array([0.1, 0.0, 0.0, 0.0])
    eps_delta: np.ndarray = np.array([65.0, 40.0, 0.0, 0.0])
    eps_inf: float = 4.0
    sigma: float = 2.0
    tau: np.ndarray = np.array([7.96e-12, 1.592e-9, 159.155e-6, 5.305e-3])


class ColeColeModel(DielectricModel):
    """Model for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: ColeColeParamters) -> None:
        self.__eps_inf = parameters.eps_inf
        self.__sigma = parameters.sigma
        self.__alpha = parameters.alpha
        self.__eps_delta = parameters.eps_delta
        self.__tau = parameters.tau

    def permitivity(self, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        divisor = 1 + (1j * omega * self.__tau) ** (1 - self.__alpha)
        eps_dispersion = self.__eps_inf + np.sum(self.__eps_delta / divisor)

        if omega == 0:
            return eps_dispersion * self.e0 + 0j

        return self.e0 * eps_dispersion + self.__sigma / (1j * omega)

    def conductivity(self, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        if omega == 0:
            return self.__sigma + 0j

        return np.conjugate(1j * omega * self.permitivity(omega=omega))


class ColeColeFourModelFactory():
    @classmethod
    def create(cls, material: Material) -> 'DielectricModel':
        material_parameters = {Material.CSF: CerebroSpinalFluidParameters,
                               Material.WHITE_MATTER: WhiteMatterParameters,
                               Material.GRAY_MATTER: GrayMatterParameters}
        return ColeColeModel(material_parameters[material])
