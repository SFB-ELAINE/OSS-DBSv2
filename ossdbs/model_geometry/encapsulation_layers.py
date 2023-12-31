from dataclasses import dataclass
from ossdbs.dielectric_model import (dielectric_models,
                                     default_dielectric_parameters,
                                     dielectric_model_parameters)


@dataclass
class EncapsulationLayer:
    """Class to store information about encapsulation layer
    """

    name: str
    material: str = None
    dielectric_model: str = None
    dielectric_parameters: dict = None
    max_h: float = 1e6

    @property
    def dielectric_properties(self):
        if self.dielectric_model is None:
            raise ValueError("You have not specified a dielectric model for the encapsulation layer {}".format(self.name))
        if self.dielectric_parameters is None:
            parameters = default_dielectric_parameters[self.dielectric_model][self.material]
        else:
            parameters = dielectric_model_parameters[self.dielectric_model](**self.dielectric_parameters)
        return dielectric_models[self.dielectric_model](parameters)
