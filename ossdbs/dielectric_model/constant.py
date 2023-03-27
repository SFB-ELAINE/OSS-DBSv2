
from dataclasses import dataclass
from ossdbs.dielectric_model import DielectricModel
from ossdbs.materials import Material


@dataclass
class ConstantDielectricModel:
    permitivity: complex
    conductivity: complex


class WhiteMatterModel(ConstantDielectricModel):
    permitivity: complex = 6.98e4
    conductivity: complex = 6.26e-2


class GrayMatterModel(ConstantDielectricModel):
    permitivity: complex = 1.64e5
    conductivity: complex = 9.88e-2


class CerebroSpinalFluidModel(ConstantDielectricModel):
    permitivity: complex = 1.09e2
    conductivity: complex = 2.0


class BloodModel(ConstantDielectricModel):
    permitivity: complex = 5.26e3
    conductivity: complex = 7e-1


class ConstantModel(DielectricModel):
    """Constant model for the dielectric spectrum of tissues."""
    MODELS = {Material.BLOOD: BloodModel,
              Material.WHITE_MATTER: WhiteMatterModel,
              Material.GRAY_MATTER: GrayMatterModel,
              Material.CSF: CerebroSpinalFluidModel
              }

    def conductivity(self, material: Material, omega: float) -> complex:
        """Return the conductivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].conductivity()

    def permitivity(self, material: Material, omega: float) -> complex:
        """Return the permitivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        return self.MODELS[material].permitivity()
