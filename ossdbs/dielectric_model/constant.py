
from dataclasses import dataclass
from ossdbs.dielectric_model import DielectricModel
from scipy.constants import epsilon_0 as e0


@dataclass
class MaterialConstantModel:
    """
    permittivity: complex
        Relative permittivity
    conductivity: complex
    """

    permittivity: complex
    conductivity: complex

    def complex_permittivity(self, omega: float) -> complex:
        """Calculate the permittivity by the angular frequency omega.

        Parameters
        ----------
        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex permittivity.
        """
        if omega == 0:
            return self.permittivity * e0 + 0j

        return e0 * self.permittivity + self.conductivity / (1j * omega)

    def complex_conductivity(self, omega: float) -> complex:
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
            return self.conductivity + 0j

        return 1j * omega * self.complex_permittivity(omega=omega)


class ConstantModel(DielectricModel):
    """Constant model for the dielectric spectrum of tissues.

    TODO: where are these values from?
    """
    MODELS = {"Blood": MaterialConstantModel(5.26e3, 7e-1),
              "White matter": MaterialConstantModel(6.98e4, 6.26e-2),
              "Gray matter": MaterialConstantModel(1.64e5, 9.88e-2),
              "CSF": MaterialConstantModel(1.09e2, 2.0),
              "Unknown": MaterialConstantModel(1.64e5, 9.88e-2),
              }
