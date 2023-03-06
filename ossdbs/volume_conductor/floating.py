
from ossdbs.contacts import Contacts
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.volume_conductor.volume_conductor_model import Solution
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


class VolumeConductorFloating(VolumeConductor):
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
        space = self.__create_space(contacts)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        sigma = self.conductivity.distribution(frequency)
        bilinear_form = self.__bilinear_form(sigma, space, contacts)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        components = solution.components[2:]
        floating_values = {contact.name: component.vec[0]
                           for contact, component
                           in zip(contacts.floating_contacts(), components)}

        potential = solution.components[0]
        current_density = sigma * ngsolve.grad(potential)

        return Solution(potential=potential,
                        current_density=current_density,
                        conductivity=sigma,
                        floating_values=floating_values,
                        frequency=frequency)

    def __create_space(self, contacts):
        boundaries = [contact.name for contact in contacts.active_contacts()]
        spaces_active = [self.mesh.h1_space(boundaries),
                         self.mesh.surfacel2_space(boundaries)]
        spaces_floating = [self.mesh.number_space()
                           for _ in contacts.floating_contacts()]
        spaces = spaces_active + spaces_floating
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    def __bilinear_form(self, sigma, space, contacts) -> ngsolve.BilinearForm:
        trial = space.TrialFunction()
        test = space.TestFunction()
        u, lam = trial[:2]
        v, mu = test[:2]
        bilinear_form = ngsolve.BilinearForm(space)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        boundaries = [contact.name for contact in contacts.floating_contacts()]
        for (ufix, vfix, boundary) in zip(trial[2:], test[2:], boundaries):
            bilinear_form += u * mu + v * lam * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form
