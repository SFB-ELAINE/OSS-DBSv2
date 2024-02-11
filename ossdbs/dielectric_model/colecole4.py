# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass

import numpy as np
from scipy.constants import epsilon_0 as e0

from .dielectric_model import DielectricModel


@dataclass
class ColeColeParameters:
    """Parameters for the 4-Cole-Cole model suggested by Gabriel et al."""

    eps_inf: float
    sigma: float
    alpha: np.ndarray  # 4 entries
    eps_delta: np.ndarray  # 4 entries
    tau: np.ndarray  # 4 entries

    def __post_init__(self):
        """Cast input to NumPy array."""
        self.alpha = np.array(self.alpha)
        self.eps_delta = np.array(self.eps_delta)
        self.tau = np.array(self.tau)

    def assert_order(self, order: int) -> bool:
        """Assert correct length of input."""
        return (
            (len(self.alpha) == order)
            and (len(self.tau) == order)
            and (len(self.eps_delta) == order)
        )


WhiteMatterColeCole4Default = ColeColeParameters(
    alpha=np.array([0.1, 0.1, 0.3, 0.02]),
    eps_delta=np.array([32.0, 100.0, 4.0e4, 3.5e7]),
    eps_inf=4.0,
    sigma=0.02,
    tau=np.array([7.958e-12, 7.958e-9, 53.052e-6, 7.958e-3]),
)

GrayMatterColeCole4Default = ColeColeParameters(
    alpha=np.array([0.1, 0.15, 0.22, 0.0]),
    eps_delta=np.array([45.0, 400.0, 2.0e5, 4.5e7]),
    eps_inf=4.0,
    sigma=0.02,
    tau=np.array([7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3]),
)

CSFColeCole4Default = ColeColeParameters(
    alpha=np.array([0.1, 0.0, 0.0, 0.0]),
    eps_delta=np.array([65.0, 40.0, 0.0, 0.0]),
    eps_inf=4.0,
    sigma=2.0,
    tau=np.array([7.96e-12, 1.592e-9, 159.155e-6, 5.305e-3]),
)


BloodColeCole4Default = ColeColeParameters(
    alpha=np.array([0.1, 0.1, 0.0, 0.0]),
    eps_delta=np.array([56.0, 5200.0, 0.0, 0.0]),
    eps_inf=4.0,
    sigma=0.7,
    tau=np.array([8.38e-12, 132.63e-9, 0, 0]),
)

default_cole_cole4_parameters = {
    "Gray matter": GrayMatterColeCole4Default,
    "Unknown": GrayMatterColeCole4Default,
    "White matter": WhiteMatterColeCole4Default,
    "CSF": CSFColeCole4Default,
    "Blood": BloodColeCole4Default,
}


class ColeCole4Model(DielectricModel):
    """Dielectric model with 4 Cole-Cole dispersions."""

    def __init__(self, parameters: ColeColeParameters):
        self._parameters = parameters
        if not self._parameters.assert_order(4):
            raise ValueError(
                "ColeCole4Model requires information about four dispersions"
            )

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
        divisor = 1 + (1j * omega * self._parameters.tau) ** (
            1 - self._parameters.alpha
        )
        eps_dispersion = self._parameters.eps_inf + np.sum(
            self._parameters.eps_delta / divisor
        )

        if omega == 0:
            return eps_dispersion * e0 + 0j

        return e0 * eps_dispersion + self._parameters.sigma / (1j * omega)

    @property
    def static_conductivity(self) -> float:
        """Return conductivity at DC / 0 Hz."""
        return self._parameters.sigma
