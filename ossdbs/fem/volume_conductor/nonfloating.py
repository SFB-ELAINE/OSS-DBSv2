from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.stimulation_signals import FrequencyDomainSignal
import ngsolve
from ossdbs.model_geometry import ModelGeometry
from ossdbs.fem.solver import Solver
from .conductivity import ConductivityCF


class VolumeConductorNonFloating(VolumeConductor):
    """Model for representing a volume conductor which evaluates the potential.
    """

    def __init__(self,
                 geometry: ModelGeometry,
                 conductivity: ConductivityCF,
                 solver: Solver,
                 order: int,
                 meshing_parameters: dict,
                 frequency_domain_signal: FrequencyDomainSignal) -> None:
        super().__init__(geometry, conductivity, solver, order, meshing_parameters, frequency_domain_signal)
        boundaries = [contact.name for contact in self.contacts.active]
        self._space = self.h1_space(boundaries=boundaries)
        self._potential = ngsolve.GridFunction(space=self._space)
        self._floating_values = {}

    def compute_solution(self,
                         frequency: float) -> ngsolve.comp.GridFunction:
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
        self._sigma = self.conductivity_cf(self.mesh, frequency)

        # update boundary condition values
        boundary_values = self.contacts.voltages
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        self._potential.Set(coefficient=coefficient,
                            VOL_or_BND=ngsolve.BND)

        u = self._space.TrialFunction()
        v = self._space.TestFunction()
        bilinear_form = ngsolve.BilinearForm(space=self._space, symmetric=True)
        bilinear_form += self._sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        linear_form = ngsolve.LinearForm(space=self._space)
        self.solver.bvp(bilinear_form, linear_form, self._potential)
