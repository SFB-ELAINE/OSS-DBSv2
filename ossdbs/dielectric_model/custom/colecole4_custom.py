
import json
import os
from ossdbs.dielectric_model import DielectricModel
import numpy as np

from ossdbs.materials import Material


class ColeColeModel():
    """Model for the dielectric spectrum of a tissue using the
    Cole-Cole equation"""

    e0 = 8.8541878128e-12

    def __init__(self,
                 eps_inf: float,
                 sigma: float,
                 alpha: np.ndarray,
                 eps_delta: np.ndarray,
                 tau: np.ndarray) -> None:
        self.eps_inf = eps_inf
        self.sigma = sigma
        self.alpha = alpha
        self.eps_delta = eps_delta
        self.tau = tau

    def permitivity(self, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        divisor = 1 + (1j * omega * self.tau) ** (1 - self.alpha)
        eps_dispersion = self.eps_inf + np.sum(self.eps_delta / divisor)

        if omega == 0:
            return eps_dispersion * self.e0 + 0j

        return self.e0 * eps_dispersion + self.sigma / (1j * omega)

    def conductivity(self, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        if omega == 0:
            return self.sigma + 0j

        return 1j * omega * self.permitivity(omega=omega)


class WhiteMatterModel(ColeColeModel):

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'WhiteMatterColeCole4Model.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.alpha = np.array(parameters['Alpha'])
        self.eps_delta = np.array(parameters['EpsilonDelta'])
        self.eps_inf = np.array(parameters['EpsilonInfinite'])
        self.sigma = np.array(parameters['Sigma'])
        self.tau = np.array(parameters['Tau'])


class GrayMatterModel(ColeColeModel):

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'GrayMatterColeCole4Model.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.alpha = np.array(parameters['Alpha'])
        self.eps_delta = np.array(parameters['EpsilonDelta'])
        self.eps_inf = np.array(parameters['EpsilonInfinite'])
        self.sigma = np.array(parameters['Sigma'])
        self.tau = np.array(parameters['Tau'])


class CerebroSpinalFluidModel(ColeColeModel):

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'CerebroSpinalFluidColeCole4Model.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.alpha = np.array(parameters['Alpha'])
        self.eps_delta = np.array(parameters['EpsilonDelta'])
        self.eps_inf = np.array(parameters['EpsilonInfinite'])
        self.sigma = np.array(parameters['Sigma'])
        self.tau = np.array(parameters['Tau'])


class BloodModel(ColeColeModel):

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'BloodColeCole4Model.json')

        with open(path, 'r') as json_file:
            parameters = json.load(json_file)

        self.alpha = np.array(parameters['Alpha'])
        self.eps_delta = np.array(parameters['EpsilonDelta'])
        self.eps_inf = np.array(parameters['EpsilonInfinite'])
        self.sigma = np.array(parameters['Sigma'])
        self.tau = np.array(parameters['Tau'])


class ColeCole4ModelCustom(DielectricModel):
    """Cole Cole model for the dielectric spectrum of tissues."""
    MODELS = {Material.BLOOD: BloodModel,
              Material.WHITE_MATTER: WhiteMatterModel,
              Material.GRAY_MATTER: GrayMatterModel,
              Material.CSF: CerebroSpinalFluidModel}

    def conductivity(self, material: Material, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].conductivity(omega)

    def permitivity(self, material: Material, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        return self.MODELS[material].permitivity(omega)
