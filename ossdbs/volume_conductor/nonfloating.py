
from ossdbs.contacts import Contacts
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.volume_conductor.volume_conductor_model import Solution
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


class VolumeConductorNonFloating(VolumeConductor):
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh
    conductivity : Conductivity
    contacts : ContactCollection
    solver : Solver
    """

    def __init__(self,
                 mesh: Mesh,
                 conductivity: Conductivity,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.solver = solver

    def compute_solution(self,
                         frequency: float,
                         contacts: Contacts) -> ngsolve.comp.GridFunction:
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
        complex_data = self.mesh.is_complex()
        sigma = self.conductivity.distribution(frequency, complex_data)
        boundaries = [contact.name for contact in contacts.active_contacts()]
        h1_space = self.mesh.h1_space(boundaries=boundaries)
        finite_elements_space = ngsolve.FESpace(spaces=[h1_space])
        space = ngsolve.CompressCompound(finite_elements_space)
        boundary_values = contacts.voltage_values()

        coefficient = self.mesh.boundary_coefficients(boundary_values)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        u = space.TrialFunction()[0]
        v = space.TestFunction()[0]
        bilinear_form = ngsolve.BilinearForm(space=space, symmetric=True)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        potential = solution.components[0]
        current_density = sigma * ngsolve.grad(potential)

        return Solution(potential=potential,
                        current_density=current_density,
                        conductivity=sigma,
                        floating_values={},
                        frequency=frequency)
