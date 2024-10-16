# Copyright 2023, 2024 Jan Philipp Payonk, Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import ngsolve

from ossdbs.fem.solver import Solver
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.model_geometry import ModelGeometry

from .conductivity import ConductivityCF


class VolumeConductorFloatingImpedance(VolumeConductor):
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
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        solution = ngsolve.GridFunction(space=self._space)
        solution.components[0].Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        self._sigma = self.conductivity_cf(self.mesh, frequency)
        bilinear_form = self.__bilinear_form(self._sigma, self._space)
        linear_form = ngsolve.LinearForm(space=self._space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        components = solution.components[1:]
        self._floating_values = {
            contact.name: component.vec[0]
            for (contact, component) in zip(self.contacts.floating, components)
        }

    def __create_space(self):
        boundaries = [contact.name for contact in self.contacts.active]

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
