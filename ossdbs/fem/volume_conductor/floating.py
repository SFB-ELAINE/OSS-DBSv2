
from ossdbs.fem.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.fem.volume_conductor.volume_conductor_model import Solution
import ngsolve


class VolumeConductorFloating(VolumeConductor):
    """Volume conductor with floating conductors
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

        boundary_values = self.model_geometry.contacts.voltage_values()
        coefficient = self.mesh.boundary_coefficients(boundary_values)
        space = self.__create_space()
        solution = ngsolve.GridFunction(space=space)
        solution.components[0].Set(coefficient=coefficient,
                                   VOL_or_BND=ngsolve.BND)
        complex_data = self.mesh.is_complex
        sigma = self.conductivity.distribution(frequency, complex_data)
        bilinear_form = self.__bilinear_form(sigma, space)
        linear_form = ngsolve.LinearForm(space=space)

        self.solver.bvp(bilinear_form, linear_form, solution)
        components = solution.components[2:]
        floating_values = {contact.name: component.vec[0]
                           for contact, component
                           in zip(self.model_geometry.contacts.floating(), components)}

        potential = solution.components[0]
        current_density = sigma * ngsolve.grad(potential)

        return Solution(potential=potential,
                        current_density=current_density,
                        conductivity=sigma,
                        floating_values=floating_values,
                        frequency=frequency)

    def __create_space(self):
        boundaries = [contact.name for contact in self.model_geometry.contacts.active()]
        spaces_active = [self.mesh.h1_space(boundaries),
                         self.mesh.surfacel2_space(boundaries)]
        spaces_floating = [self.mesh.number_space()
                           for _ in self.model_geometry.contacts.floating()]
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
        boundaries = [contact.name for contact in self.model_geometry.contacts.floating()]
        for (ufix, vfix, boundary) in zip(trial[2:], test[2:], boundaries):
            bilinear_form += (u * mu + v * lam) * ngsolve.ds(boundary)
            bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds(boundary)

        return bilinear_form
