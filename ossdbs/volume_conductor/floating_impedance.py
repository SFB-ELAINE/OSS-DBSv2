
from ossdbs.electrodes.contacts import Contacts
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.volume_conductor.volume_conductor_model import Solution
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


class VolumeConductorFloatingImpedance(VolumeConductor):
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

        boundary_values = contacts.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space()
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        complex_data = self.mesh.is_complex()
        sigma = self.conductivity.distribution(frequency, complex_data)
        bilinear_form = self.__bilinear_form(sigma, space, contacts)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        components = solution.components[1:]
        floating_values = {contact.name: component.vec[0]
                           for (contact, component)
                           in zip(contacts.floating(), components)}

        potential = solution.components[0]
        current_density = sigma * ngsolve.grad(potential)

        return Solution(potential=potential,
                        current_density=current_density,
                        conductivity=sigma,
                        floating_values=floating_values,
                        frequency=frequency)

    def __create_space(self, contacts):
        boundaries = [contact.name for contact in contacts.active()]
        h1_space = self.mesh.h1_space(boundaries=boundaries)
        number_spaces = [self.mesh.number_space()
                         for _ in contacts.floating_contacts()]
        spaces = [h1_space] + number_spaces
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    @staticmethod
    def __bilinear_form(sigma, space, contacts):
        bilinear_form = ngsolve.BilinearForm(space)
        trial = space.TrialFunction()
        test = space.TestFunction()
        u = trial[0]
        v = test[0]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        surface_impedances = contacts.floating_impedance_values()
        boundaries = [contact.name for contact in contacts.floating_contacts()]
        for (ufix, vfix, boundary) in zip(trial[1:], test[1:], boundaries):
            a = ngsolve.CoefficientFunction(1 / surface_impedances[boundary])
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)

        return bilinear_form
