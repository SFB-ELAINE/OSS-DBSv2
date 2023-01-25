
from ossdbs.boundaries import BoundaryCollection
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
from ossdbs.solver import Solver
import ngsolve


class WeakFormLaplace:

    def __init__(self, mesh: Mesh, boundaries: BoundaryCollection) -> None:
        self.__mesh = mesh
        self.__boundaries = boundaries

    def equation_terms(self, sigma: ngsolve.CoefficientFunction) -> tuple:
        active_boundaries = self.__boundaries.active_contacts()
        h1_space = self.__mesh.h1_space(boundaries=active_boundaries)
        finite_elements_space = ngsolve.FESpace(spaces=[h1_space])
        space = ngsolve.CompressCompound(finite_elements_space)
        boundary_values = self.__boundaries.voltage_values()
        coefficient = self.__mesh.boundary_coefficients(boundary_values)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        u = space.TrialFunction()[0]
        v = space.TestFunction()[0]
        bilinear_form = ngsolve.BilinearForm(space=space, symmetric=True)
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        linear_form = ngsolve.LinearForm(space=space)
        return solution, bilinear_form, linear_form


class WeakFormFloating:

    def __init__(self, mesh: Mesh, boundaries: BoundaryCollection) -> None:
        self.__mesh = mesh
        self.__boundaries = boundaries

    def equation_terms(self, sigma: ngsolve.CoefficientFunction) -> tuple:

        boundary_values = self.__boundaries.voltage_values()
        coefficient = self.__mesh.boundary_coefficients(boundary_values)
        space = self.__create_space(self.__boundaries)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)
        return solution, bilinear_form, linear_form

    def __create_space(self):
        active_contacts = self.__boundaries.active_contacts()
        spaces_active = [self.__mesh.h1_space(active_contacts),
                         self.__mesh.surfacel2_space(active_contacts)]
        spaces_floating = [self.__mesh.number_space()
                           for _ in self.__boundaries.floating_contacts()]
        spaces = spaces_active + spaces_floating
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(finite_elements_space)

    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial_functions = space.TrialFunction()
        test_functions = space.TestFunction()
        u, lam = trial_functions[:2]
        v, mu = test_functions[:2]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        boundaries = self.__boundaries.floating_contacts()
        for index, boundary in enumerate(boundaries, 2):
            ufix = trial_functions[index]
            vfix = test_functions[index]
            bilinear_form += u * mu + v * lam * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form


class WeakFormFloatingImpedance:

    def __init__(self, mesh: Mesh, boundaries: BoundaryCollection) -> None:
        self.__mesh = mesh
        self.__boundaries = boundaries

    def equation_terms(self, sigma: ngsolve.CoefficientFunction) -> tuple:
        boundary_values = self.__boundaries.voltage_values()
        coefficient = self.__mesh.boundary_coefficients(boundary_values)
        space = self.__create_space()
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)
        return solution, bilinear_form, linear_form

    def __create_space(self):
        h1_space = self.__mesh.h1_space(self.__boundaries.active_contacts())
        number_spaces = [self.mesh.number_space()
                         for _ in self.__boundaries.floating_contacts()]
        spaces = [h1_space] + number_spaces
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(finite_elements_space)

    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial_functions = space.TrialFunction()
        test_functions = space.TestFunction()
        u = trial_functions[0]
        v = test_functions[0]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        surface_impedances = self.__boundaries.floating_impedances_values()
        boundaries = self.__boundaries.floating_contacts()
        for index, boundary in enumerate(boundaries, 1):
            a = ngsolve.CoefficientFunction(1 / surface_impedances[boundary])
            ufix = trial_functions[index]
            vfix = test_functions[index]
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)

        return bilinear_form


class VolumeConductor():
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh

    conductivity : Conductivity
    """

    def __init__(self,
                 mesh: Mesh,
                 conductivity: Conductivity,
                 boundaries: BoundaryCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.boundaries = boundaries
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Parameters
        ----------
        boundaries : dict
            Dictionary containing the values for accordingly boundaries.
        frequency : float
            Frequency [Hz] of the input signal.

        Returns
        -------
        tuple
            potential and current density of volume conductor
        """

        sigma = self.conductivity.conductivity(frequency)
        active_boundaries = self.boundaries.active_contacts()
        h1_space = self.mesh.h1_space(boundaries=active_boundaries)
        finite_elements_space = ngsolve.FESpace(spaces=[h1_space])
        space = ngsolve.CompressCompound(finite_elements_space)
        boundary_values = self.boundaries.voltage_values()
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
        print(type(solution.components))
        print(solution.components[1:])
        potential = solution.components[0]

        field = ngsolve.grad(potential)
        current_density = sigma * field
        power = ngsolve.Integrate(field * ngsolve.Conj(current_density),
                                  self.mesh.ngsolvemesh())
        impedance = 1 / power

        return potential, current_density, impedance


class VolumeConductorFloating():
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh

    conductivity : Conductivity
    """

    def __init__(self,
                 mesh: Mesh,
                 conductivity: Conductivity,
                 boundaries: BoundaryCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.boundaries = boundaries
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Parameters
        ----------
        boundaries : dict
            Dictionary containing the values for accordingly boundaries.
        frequency : float
            Frequency [Hz] of the input signal.

        Returns
        -------
        tuple
            potential and current density of volume conductor
        """

        sigma = self.conductivity.conductivity(frequency)

        boundary_values = self.boundaries.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space(self.boundaries)
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)
        self.solver.bvp(bilinear_form=bilinear_form,
                        linear_form=linear_form,
                        grid_function=solution)

        potential = solution.components[0]
        floating_potentials = [comp.vec[0] for comp in solution.components[2:]]

        return potential, floating_potentials

    def __create_space(self):
        active_contacts = self.boundaries.active_contacts()
        spaces_active = [self.mesh.h1_space(active_contacts),
                         self.mesh.surfacel2_space(active_contacts)]
        spaces_floating = [self.mesh.number_space()
                           for _ in self.boundaries.floating_contacts()]
        spaces = spaces_active + spaces_floating
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(finite_elements_space)

    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial_functions = space.TrialFunction()
        test_functions = space.TestFunction()
        u, lam = trial_functions[:2]
        v, mu = test_functions[:2]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        boundaries = self.boundaries.floating_contacts()
        for index, boundary in enumerate(boundaries, 2):
            ufix = trial_functions[index]
            vfix = test_functions[index]
            bilinear_form += u * mu + v * lam * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form


class VolumeConductorFloatingImpedance():
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh

    conductivity : Conductivity
    """

    def __init__(self,
                 mesh: Mesh,
                 conductivity: Conductivity,
                 boundaries: BoundaryCollection,
                 solver: Solver) -> None:
        self.conductivity = conductivity
        self.mesh = mesh
        self.boundaries = boundaries
        self.solver = solver

    def potential(self, frequency: float) -> ngsolve.comp.GridFunction:
        """Evaluate electrical potential of volume conductor.

        Parameters
        ----------
        boundaries : dict
            Dictionary containing the values for accordingly boundaries.
        frequency : float
            Frequency [Hz] of the input signal.

        Returns
        -------
        tuple
            potential and current density of volume conductor
        """

        sigma = self.__conductivity.conductivity(frequency)

        boundary_values = self.boundaries.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space()
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form=bilinear_form,
                        linear_form=linear_form,
                        grid_function=solution)

        potential = solution.components[0]
        floating_potentials = [comp.vec[0] for comp in solution.components[1:]]

        return potential, floating_potentials

    def __create_space(self):
        h1_space = self.mesh.h1_space(self.boundaries.active_contacts())
        number_spaces = [self.mesh.number_space()
                         for _ in self.boundaries.floating_contacts()]
        spaces = [h1_space] + number_spaces
        finite_elements_space = ngsolve.FESpace(spaces=spaces)
        return ngsolve.CompressCompound(finite_elements_space)

    def __bilinear_form(self, sigma, space):
        bilinear_form = ngsolve.BilinearForm(space)
        trial_functions = space.TrialFunction()
        test_functions = space.TestFunction()
        u = trial_functions[0]
        v = test_functions[0]
        bilinear_form += sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        surface_impedances = self.boundaries.floating_impedances_values()
        boundaries = self.boundaries.floating_contacts()
        for index, boundary in enumerate(boundaries, 1):
            a = ngsolve.CoefficientFunction(1 / surface_impedances[boundary])
            ufix = trial_functions[index]
            vfix = test_functions[index]
            bilinear_form += a * (u - ufix) * (v - vfix) * ngsolve.ds(boundary)

        return bilinear_form
