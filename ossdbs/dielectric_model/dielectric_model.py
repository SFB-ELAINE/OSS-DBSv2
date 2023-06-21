import numpy as np
from abc import ABC
from ossdbs.utils.materials import MATERIALS


class DielectricModel(ABC):
    """Model for the dielectric spectrum of tissues.

    TODO explain how to add another model
    """

    @property
    def material_models(self) -> dict:
        return self.MODELS

    def permittivity(self, material: str, omega: float) -> float:
        """Calculate the permittivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        float
            Permittivity.
        """
        
        return np.real(self.MODELS[material].complex_permittivity(omega))

    def complex_permittivity(self, material: str, omega: float) -> float:
        """Calculate the permittivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex permittivity.
        """
        return self.MODELS[material].complex_permittivity(omega)

    def conductivity(self, material: str, omega: float) -> float:
        """Calculate the conductivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        float
            Conductivity.
        """
        return np.real(self.MODELS[material].complex_conductivity(omega))

    def complex_conductivity(self, material: str, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].complex_conductivity(omega)
