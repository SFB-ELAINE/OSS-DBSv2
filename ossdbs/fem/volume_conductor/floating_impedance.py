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
        raise NotImplementedError(
            "This module needs to be fixed before "
            "productive release. A proper implementation of the "
            "complete electrode model is required."
        )
        _logger.debug("Save surface impedance boundaries")
        self._surface_impedance_floating_boundaries = []
        for contact in self.contacts.floating:
            if contact.surface_impedance_model is not None:
                self._surface_impedance_floating_boundaries.append(contact.name)
            else:
                _logger.warning(
                    f"Contact {contact.name} ignored because no "
                    "surface impedance model given."
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

        # update boundary condition values
        _logger.debug("Assign potential values")
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        self._solution.components[0].Set(
            coefficient=coefficient, VOL_or_BND=ngsolve.BND
        )

        self._surface_impedances = self.contacts.get_surface_impedances(
            frequency, is_complex=self.is_complex
        )

        _logger.debug("Bilinear form")
        bilinear_form = self.__bilinear_form(self._sigma, self._space, frequency)
        _logger.debug("Linear form")
        linear_form = self.__linear_form(self._space)
        _logger.debug("Solve BVP")
        self.solver.bvp(bilinear_form, linear_form, self._solution)
        _logger.debug("Get floating values")
        self._update_floating_voltages()

    def __create_space(self) -> ngsolve.FESpace:
        h1_space = self.h1_space(boundaries=None, is_complex=self.is_complex)
        number_spaces = [
            self.number_space() for _ in self._surface_impedance_floating_boundaries
        ]
        spaces = [h1_space, *number_spaces]
        return ngsolve.FESpace(spaces=spaces)
        # very slow?!
        # return ngsolve.CompressCompound(fespace=finite_elements_space)

    def __bilinear_form(self, sigma, space, frequency) -> ngsolve.BilinearForm:
        """Bilinear form."""
        bilinear_form = ngsolve.BilinearForm(space)
        trial = space.TrialFunction()
        test = space.TestFunction()
        u = trial[0]
        v = test[0]
        _logger.debug(f"Creating space with {len(test)} subspaces")

        # sum up all potentials to set sum to zero in the end
        sum_u = None
        sum_v = None

        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx

        # add surface impedances
        for ufix, vfix, boundary in zip(
            trial[1:],
            test[1:],
            self._surface_impedance_floating_boundaries,
        ):
            print(1.0 / self._surface_impedances[boundary])
            ys = ngsolve.CoefficientFunction(1.0 / self._surface_impedances[boundary])
            bilinear_form += ys * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)
            if sum_u is None and sum_v is None:
                sum_u = ufix
                sum_v = vfix
            else:
                sum_u += ufix
                sum_v += vfix
        # enforces that sum of potentials is zero
        bilinear_form += sum_u * test[-1] * ngsolve.dx
        bilinear_form += sum_v * trial[-1] * ngsolve.dx

        return bilinear_form

    def __linear_form(self, space) -> ngsolve.LinearForm:
        v = space.TestFunction()[0]
        contacts = list(self.contacts.active)
        contacts_floating = list(self.contacts.floating)
        contacts.extend(contacts_floating)
        f = ngsolve.LinearForm(space=self._space)
        for contact in contacts:
            # account for mm as length unit
            # TODO check here
            length = ngsolve.Integrate(
                ngsolve.CoefficientFunction(1e-3) * ngsolve.ds(contact.name),
                self.mesh.ngsolvemesh,
            )
            f += contact.current / length * v * ngsolve.ds(contact.name)
        return f

    def _update_floating_voltages(self) -> None:
        """Set contact voltages with floating values."""
        components = self._solution.components[1:]
        floating_values = {
            contact: component.vec[0]
            for (contact, component) in zip(
                self._surface_impedance_floating_boundaries,
                components,
            )
        }
        for contact in self.contacts.floating:
            if contact.name in floating_values:
                contact.voltage = floating_values[contact.name]
            _logger.debug(
                f"""Contact {contact.name} updated
                    with floating potential {contact.voltage}"""
            )
