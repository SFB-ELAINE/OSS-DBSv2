
from dataclasses import dataclass
from ossdbs.dielectric_model import DielectricModel
from ossdbs.materials import Material


@dataclass
class ConstantDielectricModel:
    permitivity: complex
    conductivity: complex


class ConstantModel(DielectricModel):
    """Constant model for the dielectric spectrum of tissues."""
    MODELS = {Material.BLOOD: ConstantDielectricModel(5.26e3, 7e-1),
              Material.WHITE_MATTER: ConstantDielectricModel(6.98e4, 6.26e-2),
              Material.GRAY_MATTER: ConstantDielectricModel(1.64e5, 9.88e-2),
              Material.CSF: ConstantDielectricModel(1.09e2, 2.0),
              Material.UNKNOWN: ConstantDielectricModel(1.64e5, 9.88e-2),
              }

    def conductivity(self, material: Material, omega: float) -> complex:
        """Return the conductivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].conductivity

    def permitivity(self, material: Material, omega: float) -> complex:
        """Return the permitivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        return self.MODELS[material].permitivity
