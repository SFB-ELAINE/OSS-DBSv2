
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.fem.volume_conductor.volume_conductor_model import Solution
import ngsolve


class VolumeConductorNonFloating(VolumeConductor):
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh
    conductivity : Conductivity
    solver : Solver
    """

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
        complex_data = self.mesh.is_complex
        sigma = self.conductivity.distribution(frequency, complex_data)
        boundaries = [contact.name for contact in self.model_geometry.contacts.active()]
        space = self.mesh.h1_space(boundaries=boundaries)
        boundary_values = self.model_geometry.contacts.voltage_values()

        coefficient = self.mesh.boundary_coefficients(boundary_values)
        potential = ngsolve.GridFunction(space=space)
        potential.Set(coefficient=coefficient,
                      VOL_or_BND=ngsolve.BND)
        u = space.TrialFunction()
        v = space.TestFunction()
        bilinear_form = ngsolve.BilinearForm(space=space, symmetric=True)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, potential)
        current_density = sigma * ngsolve.grad(potential)

        return Solution(potential=potential,
                        current_density=current_density,
                        conductivity=sigma,
                        floating_values={},
                        frequency=frequency)
