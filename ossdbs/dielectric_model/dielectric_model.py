# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

from abc import ABC, abstractmethod

import numpy as np


class DielectricModel(ABC):
    """Model for the dielectric spectrum of tissues.

    Notes
    -----
    To add another model, define its complex permittivity
    and add the static conductivity.
    """

    def permittivity(self, omega: float) -> float:
        """Calculate the permittivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        float
            Permittivity.
        """
        return np.real(self.complex_permittivity(omega))

    @abstractmethod
    def complex_permittivity(self, omega: float) -> complex:
        """Calculate the permittivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        complex
            Complex permittivity.
        """
        pass

    def conductivity(self, omega: float) -> float:
        """Calculate the conductivity by the angular frequency omega.

        Parameters
        ----------
        material : str
            Corresponding material.

        omega : float
            Angular frequency [1/s].

        Returns
        -------
        float
            Conductivity.
        """
        return np.real(self.complex_conductivity(omega))

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
            return self.static_conductivity + 0j

        return 1j * omega * self.complex_permittivity(omega=omega)

    @property
    @abstractmethod
    def static_conductivity(self) -> float:
        """Return conductivity at zero Hz (static limit)."""
        pass
