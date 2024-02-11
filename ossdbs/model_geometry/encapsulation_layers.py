# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass

from ossdbs.dielectric_model import (
    default_dielectric_parameters,
    dielectric_model_parameters,
    dielectric_models,
)


@dataclass
class EncapsulationLayer:
    """Class to store information about encapsulation layer."""

    name: str
    material: str = None
    dielectric_model: str = None
    dielectric_parameters: dict = None
    max_h: float = 1e6

    @property
    def dielectric_properties(self):
        """Get dielectric properties of encapsulation layer."""
        if self.dielectric_model is None:
            raise ValueError(
                f"""You have not specified a dielectric model
                for the encapsulation layer {self.name}"""
            )
        if self.dielectric_parameters is None:
            parameters = default_dielectric_parameters[self.dielectric_model][
                self.material
            ]
        else:
            parameters = dielectric_model_parameters[self.dielectric_model](
                **self.dielectric_parameters
            )
        return dielectric_models[self.dielectric_model](parameters)
