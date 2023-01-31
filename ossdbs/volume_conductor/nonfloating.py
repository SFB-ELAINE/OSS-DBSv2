
from ossdbs.volume_conductor.volume_conductor_model import VolumeConductor
from ossdbs.volume_conductor.volume_conductor_model import Potential
from ossdbs.electrode_contacts import ContactCollection
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
