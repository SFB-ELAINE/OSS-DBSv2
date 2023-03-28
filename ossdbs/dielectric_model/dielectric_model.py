
from abc import ABC, abstractmethod

from ossdbs.materials import Material


class DielectricModel(ABC):
    """Model for the dielectric spectrum of tissues."""

    @abstractmethod
    def permitivity(self, material: Material, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        pass

    @abstractmethod
    def conductivity(self, material: Material, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        pass
