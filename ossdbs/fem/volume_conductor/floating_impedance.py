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
        # Include both active contacts (Dirichlet BC) and floating+impedance contacts
        active_boundaries = [contact.name for contact in self.contacts.active]
        boundaries = active_boundaries + self._surface_impedance_floating_boundaries
        h1_space = self.h1_space(boundaries=boundaries, is_complex=self.is_complex)
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
        _logger.debug("Creating space with {len(u)} subspaces")

        # sum up all potentials to set sum to zero in the end
        sum_u = None
        sum_v = None

        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx

        # add surface impedances
        for ufix, vfix, boundary in zip(
            trial[1:],
            test[1:],
            self._surface_impedance_floating_boundaries,
            strict=True,
        ):
            a = ngsolve.CoefficientFunction(1.0 / self._surface_impedances[boundary])
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)
            if sum_u is None and sum_v is None:
                sum_u = ufix
                sum_v = vfix
            else:
                sum_u += ufix
                sum_v += vfix

        # enforces that sum of potentials is zero (only needed if no active
        # contacts). With Dirichlet BC the system is already well-posed.
        if len(self.contacts.active) == 0 and sum_u is not None:
            bilinear_form += sum_u * test[-1] * ngsolve.dx
            bilinear_form += sum_v * trial[-1] * ngsolve.dx

        return bilinear_form

    def __linear_form(self, space) -> ngsolve.LinearForm:
        test = space.TestFunction()
        f = ngsolve.LinearForm(space=self._space)

        # For active contacts, apply current to H1 space
        v_h1 = test[0]
        for contact in self.contacts.active:
            length = ngsolve.Integrate(
                ngsolve.CoefficientFunction(1e-3) * ngsolve.ds(contact.name),
                self.mesh.ngsolvemesh,
            )
            f += contact.current / length * v_h1 * ngsolve.ds(contact.name)

        # For floating+impedance contacts, apply current to NumberSpace DOF
        for idx, boundary in enumerate(self._surface_impedance_floating_boundaries):
            contact = None
            for c in self.contacts.floating:
                if c.name == boundary:
                    contact = c
                    break
            if contact is None:
                continue
            v_float = test[idx + 1]  # NumberSpace component
            length = ngsolve.Integrate(
                ngsolve.CoefficientFunction(1e-3) * ngsolve.ds(contact.name),
                self.mesh.ngsolvemesh,
            )
            f += contact.current / length * v_float * ngsolve.ds(contact.name)

        return f

    def _update_floating_voltages(self) -> None:
        """Set contact voltages with floating values."""
        components = self._solution.components[1:]
        floating_values = {
            contact: component.vec[0]
            for (contact, component) in zip(
                self._surface_impedance_floating_boundaries,
                components,
                strict=True,
            )
        }
        for contact in self.contacts.floating:
            if contact.name in floating_values:
                contact.voltage = floating_values[contact.name]
            _logger.debug(
                f"""Contact {contact.name} updated
                    with floating potential {contact.voltage}"""
            )
