
from ossdbs.dielectric_model import DielectricModel
from ossdbs.materials import Material
import numpy as np


class MaterialColeColeModel():
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

        Parameters
        ----------
        omega : float
            Angular frequency [1/s].

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

        Parameters
        ----------
        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex conductivity.
        """
        if omega == 0:
            return self.sigma + 0j

        return 1j * omega * self.permitivity(omega=omega)


WhiteMatterModel = MaterialColeColeModel(
                    alpha=np.array([0.1, 0.1, 0.3, 0.02]),
                    eps_delta=np.array([32.0, 100.0, 4.0e4, 3.5e7]),
                    eps_inf=4.0,
                    sigma=0.02,
                    tau=np.array([7.958e-12, 7.958e-9, 53.052e-6, 7.958e-3])
                    )

GrayMatterModel = MaterialColeColeModel(
                    alpha=np.array([0.1, 0.15, 0.22, 0.0]),
                    eps_delta=np.array([45.0, 400.0, 2.0e5, 4.5e7]),
                    eps_inf=4.0,
                    sigma=0.02,
                    tau=np.array([7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3])
                    )

CSFModel = MaterialColeColeModel(
                    alpha=np.array([0.1, 0.0, 0.0, 0.0]),
                    eps_delta=np.array([65.0, 40.0, 0.0, 0.0]),
                    eps_inf=4.0,
                    sigma=2.0,
                    tau=np.array([7.96e-12, 1.592e-9, 159.155e-6, 5.305e-3]),
                    )


BloodModel = MaterialColeColeModel(
                                alpha=np.array([0.1, 0.1, 0.0, 0.0]),
                                eps_delta=np.array([56.0, 5200.0, 0.0, 0.0]),
                                eps_inf=4.0,
                                sigma=0.7,
                                tau=np.array([8.38e-12, 132.63e-9, 0, 0])
                                )


class ColeCole4Model(DielectricModel):
    """Cole Cole model for the dielectric spectrum of tissues."""
    MODELS = {Material.BLOOD: BloodModel,
              Material.WHITE_MATTER: WhiteMatterModel,
              Material.GRAY_MATTER: GrayMatterModel,
              Material.CSF: CSFModel,
              Material.UNKNOWN: GrayMatterModel
              }

    def conductivity(self, material: Material, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Parameters
        ----------
        material : Material
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].conductivity(omega)

    def permitivity(self, material: Material, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Parameters
        ----------
        material : Material
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex permitivity.
        """
        return self.MODELS[material].permitivity(omega)
