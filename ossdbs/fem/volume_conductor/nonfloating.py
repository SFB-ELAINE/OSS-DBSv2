# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import ngsolve

from ossdbs.fem.solver import Solver
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.model_geometry import ModelGeometry

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductorNonFloating(VolumeConductor):
    """Model for representing a volume conductor which evaluates the potential."""

    def __init__(
        self,
        geometry: ModelGeometry,
        conductivity: ConductivityCF,
        solver: Solver,
        order: int,
        meshing_parameters: dict,
        output_path: str = "Results",
    ) -> None:
        super().__init__(
            geometry, conductivity, solver, order, meshing_parameters, output_path
        )
        boundaries = [contact.name for contact in self.contacts.active]
        self._space = self.h1_space(boundaries=boundaries, is_complex=self.is_complex)
        self._floating_values = {}
        self.update_space()

    def update_space(self):
        """Update space (e.g., if mesh changes)."""
        # only GridFunction needs to be updated
        self._potential = ngsolve.GridFunction(space=self._space)

    def compute_solution(self, frequency: float) -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Parameters
        ----------
        frequency : float
            Frequency [Hz] of the input signal.

        Returns
        -------
        Potential
            Data object representing the potential of volume conductor and
            floating values of floating contacts.
        """
        self._frequency = frequency
        _logger.debug("Get conductivity at frequency")
        self._sigma = self.conductivity_cf(self.mesh, frequency)
        _logger.debug(f"Sigma: {self._sigma}")

        # update boundary condition values
        _logger.debug("Assign potential values")
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        self._potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)

        _logger.debug("Prepare weak form")
        u = self._space.TrialFunction()
        v = self._space.TestFunction()
        _logger.debug("Bilinear form")
        # TODO symmetric even if anisotropy?
        bilinear_form = ngsolve.BilinearForm(space=self._space, symmetric=True)
        _logger.debug("Bilinear form, formulation")
        bilinear_form += self._sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx

        # TODO add surface impedance
        _logger.debug("Linear form")
        linear_form = ngsolve.LinearForm(space=self._space)
        _logger.debug("Solve BVP")
        self.solver.bvp(bilinear_form, linear_form, self._potential)
