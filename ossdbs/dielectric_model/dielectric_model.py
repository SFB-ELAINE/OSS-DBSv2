
from abc import ABC, abstractmethod


class DielectricModel(ABC):
    """Model for the dielectric spectrum of a tissue."""

    @abstractmethod
    def permitivity(self, omega: float) -> complex:
        """Calculate the permitivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """

        pass

    @abstractmethod
    def conductivity(self, omega: float) -> complex:
        """Calculate the conductivity by the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """

        pass
