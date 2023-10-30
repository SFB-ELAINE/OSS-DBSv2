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
        self._floating_values = {}
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
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        self._potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        _logger.debug("Bilinear form")
        bilinear_form = self.__bilinear_form(self._sigma, self._space)
        _logger.debug("Linear form")
        linear_form = self.__linear_form(self._space)
        _logger.debug("Solve BVP")
        self.solver.bvp(bilinear_form, linear_form, self._potential)
        _logger.debug("Get floating values")
        self._floating_values = {
                # TODO
        }

    def __create_space(self) -> ngsolve.H1 :
        boundaries = [contact.name for contact in self.contacts.active]
        if len(boundaries) == 0:
            raise RuntimeError("At least one boundary with a fixed voltage has to be specified!")
        boundaries_floating = [contact.name for contact in self.contacts.floating]
        space_field = self.h1_space(boundaries)
        plateaus = []
        for boundary in boundaries_floating:
            plateaus.append(self.mesh.ngsolvemesh.Boundaries(boundary))
        _logger.debug("Create finite element space")
        finite_element_space = ngsolve.PlateauFESpace(space_field, plateaus)
        return finite_element_space

    def __bilinear_form(self, sigma, space) -> ngsolve.BilinearForm:
        u = space.TrialFunction()
        v = space.TestFunction()
        bilinear_form = ngsolve.BilinearForm(space)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        return bilinear_form

    def __linear_form(self, space) -> ngsolve.LinearForm:
        v = space.TestFunction()
        contacts = [contact for contact in self.contacts.active]
        contacts_floating = [contact for contact in self.contacts.floating]
        contacts.extend(contacts_floating)
        f = ngsolve.LinearForm(space=self._space)
        for contact in contacts:
            length = ngsolve.Integrate(ngsolve.CoefficientFunction(1.0) * ngsolve.ds(contact.name), self.mesh.ngsolvemesh)
            f +=  contact.current / length * v * ngsolve.ds(contact.name)
        return f 
