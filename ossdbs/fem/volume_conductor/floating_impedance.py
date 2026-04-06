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
        # Only active contacts get Dirichlet BC; floating+impedance contacts
        # are coupled via Robin BC and must NOT be Dirichlet
        active_boundaries = [contact.name for contact in self.contacts.active]
        h1_space = self.h1_space(
            boundaries=active_boundaries, is_complex=self.is_complex
        )
        number_spaces = [
            self.number_space() for _ in self._surface_impedance_floating_boundaries
        ]
        spaces = [h1_space, *number_spaces]
        # If no active contacts, add a Lagrange multiplier NumberSpace
        # to enforce sum of floating potentials = 0
        self._needs_lagrange = len(self.contacts.active) == 0 and len(number_spaces) > 0
        if self._needs_lagrange:
            spaces.append(self.number_space())
        return ngsolve.FESpace(spaces=spaces)

    def __bilinear_form(self, sigma, space, frequency) -> ngsolve.BilinearForm:
        """Bilinear form."""
        bilinear_form = ngsolve.BilinearForm(space)
        trial = space.TrialFunction()
        test = space.TestFunction()
        u = trial[0]
        v = test[0]
        _logger.debug("Creating space with {len(u)} subspaces")

        n_float = len(self._surface_impedance_floating_boundaries)

        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx

        # add surface impedances
        sum_u = None
        sum_v = None
        for ufix, vfix, boundary in zip(
            trial[1 : 1 + n_float],
            test[1 : 1 + n_float],
            self._surface_impedance_floating_boundaries,
            strict=True,
        ):
            a = ngsolve.CoefficientFunction(1.0 / self._surface_impedances[boundary])
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)
            if sum_u is None:
                sum_u = ufix
                sum_v = vfix
            else:
                sum_u += ufix
                sum_v += vfix

        # Lagrange multiplier to enforce sum of floating potentials = 0
        if self._needs_lagrange and sum_u is not None:
            lam = trial[1 + n_float]
            mu = test[1 + n_float]
            bnd = self._surface_impedance_floating_boundaries[0]
            bilinear_form += sum_u * mu * ngsolve.ds(bnd)
            bilinear_form += sum_v * lam * ngsolve.ds(bnd)
            # Small regularization so the Lagrange multiplier block
            # has a non-zero diagonal (needed for Jacobi-type
            # preconditioners in iterative solvers like GMRES).
            eps = 1e-12
            bilinear_form += eps * lam * mu * ngsolve.ds(bnd)

        return bilinear_form

    def __linear_form(self, space) -> ngsolve.LinearForm:
        test = space.TestFunction()
        f = ngsolve.LinearForm(space=self._space)

        # For active contacts, apply current to H1 space
        v_h1 = test[0]
        for contact in self.contacts.active:
            length = ngsolve.Integrate(
                ngsolve.CoefficientFunction(1.0) * ngsolve.ds(contact.name),
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
                ngsolve.CoefficientFunction(1.0) * ngsolve.ds(contact.name),
                self.mesh.ngsolvemesh,
            )
            f += contact.current / length * v_float * ngsolve.ds(contact.name)

        return f

    def _update_floating_voltages(self) -> None:
        """Set contact voltages with floating values."""
        n_float = len(self._surface_impedance_floating_boundaries)
        components = self._solution.components[1 : 1 + n_float]
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
