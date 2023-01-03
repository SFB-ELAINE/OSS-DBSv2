from abc import ABC
import numpy as np


class DielectricParamters(ABC):
    alpha: np.ndarray
    eps_delta: np.ndarray
    eps_inf: float
    sigma: float
    tau: np.ndarray


class DielectricModel():
    """Model for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self, parameters: DielectricParamters) -> None:
        self.__eps_inf = parameters.eps_inf
        self.__sigma = parameters.sigma
        self.__alpha = parameters.alpha
        self.__eps_delta = parameters.eps_delta
        self.__tau = parameters.tau

    def permitivity(self, omega: float) -> complex:
        divisor = 1 + (1j * omega * self.__tau) ** (1 - self.__alpha)
        eps_dispersion = self.__eps_inf + np.sum(self.__eps_delta / divisor)

        if omega == 0:
            return eps_dispersion * self.e0 + 0j

        return self.e0 * eps_dispersion + self.__sigma / (1j * omega)

    def conductivity(self, omega: float) -> complex:

        if omega == 0:
            return self.__sigma + 0j

        return np.conjugate(1j * omega * self.permitivity(omega=omega))
