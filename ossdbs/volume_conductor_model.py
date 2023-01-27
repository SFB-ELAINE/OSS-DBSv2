
from abc import ABC, abstractmethod
from dataclasses import dataclass
from ossdbs.electrode_contacts import ContactCollection
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


@dataclass
class Potential:
    gridfunction: ngsolve.GridFunction
    frequency: float
    floating_values: dict


class VolumeConductor(ABC):

    @abstractmethod
    def potential(self, frequency: float) -> Potential:
        pass


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
                 contacts: ContactCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.contacts = contacts
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
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

        sigma = self.conductivity.distribution(frequency)
        active_boundaries = self.contacts.active()
        h1_space = self.mesh.h1_space(boundaries=active_boundaries)
        finite_elements_space = ngsolve.FESpace(spaces=[h1_space])
        space = ngsolve.CompressCompound(finite_elements_space)
        boundary_values = self.contacts.voltage_values()
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

        field = ngsolve.grad(potential)
        current_density = sigma * field
        power = ngsolve.Integrate(field * ngsolve.Conj(current_density),
                                  self.mesh.ngsolvemesh())
        impedance = 1 / power
        print(impedance)

        return Potential(gridfunction=solution.components[0],
                         floating_values={},
                         frequency=frequency)


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
                 contacts: ContactCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.contacts = contacts
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
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

        boundary_values = self.contacts.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space(self.contacts)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        sigma = self.conductivity.distribution(frequency)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        floating_results = zip(self.contacts.floating(),
                               solution.components[2:])
        floating_values = {contact: component.vec[0]
                           for (contact, component) in floating_results}

        return Potential(gridfunction=solution.components[0],
                         floating_values=floating_values,
                         frequency=frequency)

    def __create_space(self):
        active_contacts = self.contacts.active()
        spaces_active = [self.mesh.h1_space(active_contacts),
                         self.mesh.surfacel2_space(active_contacts)]
        spaces_floating = [self.mesh.number_space()
                           for _ in self.contacts.floating()]
        spaces = spaces_active + spaces_floating
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    def __bilinear_form(self, sigma, space) -> ngsolve.BilinearForm:
        trial = space.TrialFunction()
        test = space.TestFunction()
        u, lam = trial[:2]
        v, mu = test[:2]
        bilinear_form = ngsolve.BilinearForm(space)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        boundaries = self.contacts.floating()
        for (ufix, vfix, boundary) in zip(trial[2:], test[2:], boundaries):
            bilinear_form += u * mu + v * lam * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form


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
                 contacts: ContactCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.contacts = contacts
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
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

        boundary_values = self.contacts.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space()
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        sigma = self.conductivity.distribution(frequency)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        floating_results = zip(self.contacts.floating(),
                               solution.components[1:])
        floating_values = {contact: component.vec[0]
                           for (contact, component) in floating_results}

        return Potential(gridfunction=solution.components[0],
                         floating_values=floating_values,
                         frequency=frequency)

    def __create_space(self):
        h1_space = self.mesh.h1_space(self.contacts.active())
        number_spaces = [self.mesh.number_space()
                         for _ in self.boundaries.floating_contacts()]
        spaces = [h1_space] + number_spaces
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(fespace=finite_elements_space)

    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial = space.TrialFunction()
        test = space.TestFunction()
        u = trial[0]
        v = test[0]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        surface_impedances = self.contacts.floating_impedance_values()
        boundaries = self.contacts.floating()
        for (ufix, vfix, boundary) in zip(trial[1:], test[1:], boundaries):
            a = ngsolve.CoefficientFunction(1 / surface_impedances[boundary])
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)

        return bilinear_form
