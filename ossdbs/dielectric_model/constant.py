# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass

from scipy.constants import epsilon_0 as e0

from .dielectric_model import DielectricModel


@dataclass
class ConstantParameters:
    """Parameters for the constant model.
    The default values are taken from the model by Gabriel et al.
    at a frequency of 10 kHz.
    """

    permittivity: float
    conductivity: float


BloodConstantDefault = ConstantParameters(permittivity=5.25e3, conductivity=7e-1)
WhiteMatterConstantDefault = ConstantParameters(
    permittivity=1.25e4, conductivity=6.95e-2
)
GrayMatterConstantDefault = ConstantParameters(
    permittivity=2.22e4, conductivity=1.15e-1
)
CSFConstantDefault = ConstantParameters(permittivity=1.09e2, conductivity=2.0)


default_constant_parameters = {
    "Gray matter": GrayMatterConstantDefault,
    "Unknown": GrayMatterConstantDefault,
    "White matter": WhiteMatterConstantDefault,
    "CSF": CSFConstantDefault,
    "Blood": BloodConstantDefault,
}


class ConstantModel(DielectricModel):
    """Model for frequency-independent properties."""

    def __init__(self, parameters: ConstantParameters):
        self._parameters = parameters

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
            return self._parameters.permittivity * e0 + 0j

        return e0 * self._parameters.permittivity + self._parameters.conductivity / (
            1j * omega
        )

    @property
    def static_conductivity(self) -> float:
        """Returns the static conductivity."""
        return self._parameters.conductivity
