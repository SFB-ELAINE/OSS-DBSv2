import logging

import ngsolve
from ossdbs.fem.solver import Solver
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.model_geometry import ModelGeometry
from ossdbs.stimulation_signals import FrequencyDomainSignal

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductorFloating(VolumeConductor):
    """Volume conductor with floating conductors."""

    def __init__(
        self,
        geometry: ModelGeometry,
        conductivity: ConductivityCF,
        solver: Solver,
        order: int,
        meshing_parameters: dict,
        frequency_domain_signal: FrequencyDomainSignal,
    ) -> None:
        super().__init__(
            geometry,
            conductivity,
            solver,
            order,
            meshing_parameters,
            frequency_domain_signal,
        )
        _logger.debug("Create space")
        self._space = self.__create_space()
        print("Ndofs compressed fes", self._space.ndof)
        self._floating_values = {}
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
        self._potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        _logger.debug("Bilinear form")
        bilinear_form = self.__bilinear_form(self._sigma, self._space)
        _logger.debug("Linear form")
        linear_form = ngsolve.LinearForm(space=self._space)
        _logger.debug("Solve BVP")
        self.solver.bvp(bilinear_form, linear_form, self._solution)
        _logger.debug("Get floating values")
        components = self._solution.components[2:]
        self._floating_values = {
            contact.name: component.vec[0]
            for contact, component in zip(
                self.contacts.floating, components
            )
        }

    def __create_space(self):
        boundaries = [contact.name for contact in self.contacts.active]
        boundaries_nonfloating = [
            contact.name for contact in self.contacts.unused
        ]
        boundaries_nonfloating.extend(boundaries)
        spaces_field = [
            self.h1_space(boundaries),
            self.surfacel2_space(boundaries_nonfloating),
        ]
        spaces_floating = [
            self.number_space() for _ in self.contacts.floating
        ]
        spaces = spaces_field + spaces_floating
        _logger.debug("Create finite element space")
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        _logger.debug("CompressComound finite element space")
        print("Ndofs finite_element_space", finite_elements_space.ndof)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    def __bilinear_form(self, sigma, space) -> ngsolve.BilinearForm:
        trial = space.TrialFunction()
        test = space.TestFunction()
        u, lam = trial[:2]
        v, mu = test[:2]
        bilinear_form = ngsolve.BilinearForm(space)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        boundaries = [contact.name for contact in self.contacts.floating]
        for ufix, vfix, boundary in zip(trial[2:], test[2:], boundaries):
            bilinear_form += (u * mu + v * lam) * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form
