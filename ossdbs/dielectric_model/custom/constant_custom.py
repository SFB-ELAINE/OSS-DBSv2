
import json
import os
from ossdbs.dielectric_model import DielectricModel
from ossdbs.materials import Material


class WhiteMatterModel():

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'WhiteMatterConstantModel.json')

        with open(path, 'r') as json_file:
            par = json.load(json_file)

        self.permitivity = par['PermitivityReal'] + par['PermitivityImag']
        self.conductivity = par['ConductivityReal'] + par['ConductivityImag']


class GrayMatterModel():

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'GrayMatterConstantModel.json')

        with open(path, 'r') as json_file:
            par = json.load(json_file)

        self.permitivity = par['PermitivityReal'] + par['PermitivityImag']
        self.conductivity = par['ConductivityReal'] + par['ConductivityImag']


class CerebroSpinalFluidModel():

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'CerebroSpinalFluidConstantModel.json')

        with open(path, 'r') as json_file:
            par = json.load(json_file)

        self.permitivity = par['PermitivityReal'] + par['PermitivityImag']
        self.conductivity = par['ConductivityReal'] + par['ConductivityImag']


class BloodModel():

    def __init__(self) -> None:

        dir_name = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(dir_name, 'BloodConstantModel.json')

        with open(path, 'r') as json_file:
            par = json.load(json_file)

        self.permitivity = par['PermitivityReal'] + par['PermitivityImag']
        self.conductivity = par['ConductivityReal'] + par['ConductivityImag']


class ConstantModelCustom(DielectricModel):
    """Constant model for the dielectric spectrum of tissues."""
    MODELS = {Material.BLOOD: BloodModel,
              Material.WHITE_MATTER: WhiteMatterModel,
              Material.GRAY_MATTER: GrayMatterModel,
              Material.CSF: CerebroSpinalFluidModel}

    def conductivity(self, material: Material, omega: float) -> complex:
        """Return the conductivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex conductivity.
        """
        return self.MODELS[material].conductivity(omega)

    def permitivity(self, material: Material, omega: float) -> complex:
        """Return the permitivity independent of the angular frequency omega.

        Returns
        -------
        complex
            Complex permitivity.
        """
        return self.MODELS[material].permitivity(omega)
