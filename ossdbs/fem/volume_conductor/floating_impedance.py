# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import ngsolve

from ossdbs.fem.solver import Solver
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.model_geometry import ModelGeometry

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductorFloatingImpedance(VolumeConductor):
    """Volume conductor with floating conductors and surface impedances."""

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
        _logger.debug("Create space")
        self._floating_values = {}
        self.update_space()

    def update_space(self):
        """Update space (e.g., if mesh changes)."""
        self._space = self.__create_space()
        self._solution = ngsolve.GridFunction(space=self._space)
        self._potential = self._solution.components[0]

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
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        self._solution.components[0].Set(
            coefficient=coefficient, VOL_or_BND=ngsolve.BND
        )
        _logger.debug("Bilinear form")
        bilinear_form = self.__bilinear_form(self._sigma, self._space)
        _logger.debug("Linear form")
        # TODO how to impose current?
        linear_form = ngsolve.LinearForm(space=self._space)
        _logger.debug("Solve BVP")
        self.solver.bvp(bilinear_form, linear_form, self._solution)
        _logger.debug("Get floating values")
        self._update_floating_voltages()

    def __create_space(self) -> ngsolve.FESpace:
        boundaries = [contact.name for contact in self.contacts.active]
        # TODO insert condition similar to EIT?
        h1_space = self.h1_space(boundaries=boundaries)
        number_spaces = [self.number_space() for _ in self.contacts.floating]
        spaces = [h1_space, *number_spaces]
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    @staticmethod
    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial = space.TrialFunction()
        test = space.TestFunction()
        u = trial[0]
        v = test[0]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        surface_impedances = self.contacts.floating_impedance_values()
        boundaries = [contact.name for contact in self.contacts.floating]
        for ufix, vfix, boundary in zip(trial[1:], test[1:], boundaries):
            a = ngsolve.CoefficientFunction(1 / surface_impedances[boundary])
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)

        return bilinear_form

    def _update_floating_voltages(self) -> None:
        """Set contact voltages with floating values."""
        components = self._solution.components[1:]
        floating_values = {
            contact.name: component.vec[0]
            for (contact, component) in zip(self.contacts.floating, components)
        }
        for contact in self.contacts.floating:
            if contact.name in floating_values:
                contact.voltage = floating_values[contact.name]
            _logger.debug(
                f"""Contact {contact.name} updated
                    with floating potential {contact.voltage}"""
            )
