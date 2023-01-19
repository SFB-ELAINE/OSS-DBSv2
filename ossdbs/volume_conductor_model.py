
from ossdbs.conductivity import Conductivity
from ossdbs.mesh import Mesh
import ngsolve
import numpy as np


class Potential():
    """Represents the electrical potential.

    Parameters
    ----------
    gridfunction : ngsolve.GridFunction
    mesh : Mesh
    """

    def __init__(self, gridfunction: ngsolve.GridFunction, mesh: Mesh) -> None:
        self.gridfunction = gridfunction
        self.mesh = mesh

    def __add__(self, other) -> 'Potential':
        self.gridfunction.vec.data += other._gridfunction.vec.data

    def __mul__(self, value) -> 'Potential':
        self.gridfunction.vec.data *= value


class VolumeConductor():
    """Model for representing a volume conductor which evaluates the potential.

    Parameters
    ----------
    mesh : Mesh

    conductivity : Conductivity
    """

    def __init__(self, mesh: Mesh, conductivity: Conductivity) -> None:
        self.__conductivity = conductivity
        self.mesh = mesh

    # Julius: not sure if an empty dict for preconditioner_options is best idea
    # Julius: would refactor: add solver, preconditioner, etc. to __init__ of volume_conductor_model
    def potential(self, boundaries: dict, frequency: float, solver: str = "CG",
                  preconditioner: str = "bddc", preconditioner_options: dict = {"coarsetype": "local"}) \
            -> ngsolve.comp.GridFunction:
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

        conductivities = self.__conductivity.conductivity(frequency)
        # convert conductivity [S/m] to [S/mm] since mesh dimension is in mm.
        values = conductivities.data * 1e-3
        start, end = conductivities.start, conductivities.end
        if not self.mesh.is_complex():
            values = np.real(values)
        sigma = ngsolve.VoxelCoefficient(start, end, values, linear=False)
        # TODO add surface_impedance through input dictionary
        # it is a (complex) value defined on each surface
        surface_impedance = None
        # WIP: add floating electrodes
        # Missing:
        # Make sure that there is only floating or floating_impedance
        # Some calls here should be wrapped in external functions for the sake
        # of readability
        floating_electrodes = self.mesh.get_floating_electrodes()
        floating_impedance_electrodes = self.mesh.get_floating_impedance_electrodes()
        if len(floating_electrodes) > 0:
            if len(floating_impedance_electrodes) > 0:
                raise RuntimeError("You did specify both floating and floating impedance electrodes!")
            h1_space = self.mesh.h1_space()
            surfacel2_space = self.mesh.surfacel2_space()
            number_spaces = []
            for _ in range(len(floating_electrodes)):
                number_spaces.append(self.mesh.number_space())
            space = ngsolve.FESpace(h1_space, surfacel2_space, *number_spaces)
            space = ngsolve.CompressCompound(space)
            solution = ngsolve.GridFunction(space=space)
            potential = solution.components[0]
            self.bilinear_form = ngsolve.BilinearForm(space)
            self.linear_form = ngsolve.LinearForm(space)
            trial_functions = space.TrialFunction()
            test_functions = space.TestFunction()
            self.apply_weak_form_laplace_equation(trial_functions[0],
                                                  test_functions[0],
                                                  sigma)
            for idx in range(len(floating_electrodes)):
                self.apply_weak_form_floating_electrode(trial_functions[0],
                                                        trial_functions[1],
                                                        trial_functions[2 + idx],
                                                        test_functions[0],
                                                        test_functions[1],
                                                        test_functions[2 + idx],
                                                        idx)

        elif len(floating_impedance_electrodes) > 0:
            raise NotImplementedError("Surface impedance not yet implemented")
            h1_space = self.mesh.h1_space()
            number_spaces = []
            for _ in range(len(floating_impedance_electrodes)):
                number_spaces.append(self.mesh.number_space())
            space = ngsolve.FESpace(h1_space, *number_spaces)
            solution = ngsolve.GridFunction(space=space)
            potential = solution.components[0]
            self.bilinear_form = ngsolve.BilinearForm(space)
            self.linear_form = ngsolve.LinearForm(space)
            trial_functions = space.TrialFunction()
            test_functions = space.TestFunction()
            for idx in range(len(floating_impedance_electrodes)):
                self.apply_weak_form_floating_impedance_electrode(trial_functions[0],
                                                                  trial_functions[1 + idx],
                                                                  test_functions[0],
                                                                  test_functions[1 + idx],
                                                                  idx,
                                                                  surface_impedance
                                                                  )
        else:
            space = self.mesh.h1_space()
            potential = ngsolve.GridFunction(space=space)
            # ugly, I am sorry, please think about a better solution
            solution = potential
            self.bilinear_form = ngsolve.BilinearForm(space=space, symmetric=True)
            self.linear_form = ngsolve.LinearForm(space=space)
            self.apply_weak_form_laplace_equation(space.TrialFunction(),
                                                  space.TestFunction(),
                                                  sigma)

        coefficient = self.mesh.boundary_coefficients(boundaries=boundaries)

        fem_system = FEMSolver(self.bilinear_form, self.linear_form, solver, preconditioner, preconditioner_options)
        potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)
        solution.vec.data = fem_system.solve_bvp(input=solution)

        field = ngsolve.grad(potential)
        current_density = sigma * field
        power = ngsolve.Integrate(field * ngsolve.Conj(current_density),
                                  self.mesh.ngsolvemesh())
        impedance = 1 / power

        floating_potentials = []
        for idx in range(len(floating_electrodes)):
            floating_potentials.append(solution.components[2 + idx].vec[0])
        for idx in range(len(floating_impedance_electrodes)):
            floating_potentials.append(solution.components[1 + idx].vec[0])

        return potential, current_density, impedance, floating_potentials

    def apply_weak_form_floating_electrode(self,
                                           u: ngsolve.comp.FESpace.TrialFunction,
                                           lam: ngsolve.comp.FESpace.TrialFunction,
                                           ufix: ngsolve.comp.FESpace.TrialFunction,
                                           v: ngsolve.comp.FESpace.TestFunction,
                                           mu: ngsolve.comp.FESpace.TestFunction,
                                           vfix: ngsolve.comp.FESpace.TestFunction,
                                           idx: int
                                           ) -> None:

        self.bilinear_form += (u * mu + v * lam) * ngsolve.ds("floating_{}".format(idx))
        self.bilinear_form += -(ufix * mu + vfix * lam) * ngsolve.ds("floating_{}".format(idx))

        # TODO for linear form: current can be imposed

    def apply_weak_form_floating_impedance_electrode(self,
                                                     u: ngsolve.comp.FESpace.TrialFunction,
                                                     ufix: ngsolve.comp.FESpace.TrialFunction,
                                                     v: ngsolve.comp.FESpace.TestFunction,
                                                     vfix: ngsolve.comp.FESpace.TestFunction,
                                                     idx: int,
                                                     surface_impedance: ngsolve.CoefficientFunction
                                                     ) -> None:

        self.bilinear_form += 1. / surface_impedance * (u - ufix) * (v - vfix) * ngsolve.ds("floating_impedance_{}".format(idx))

        # TODO for linear form: current can be imposed

    def apply_weak_form_laplace_equation(self,
                                         u: ngsolve.comp.FESpace.TrialFunction,
                                         v: ngsolve.comp.FESpace.TestFunction,
                                         sigma: ngsolve.fem.CoefficientFunction
                                         ) -> None:

        equation = sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx
        self.bilinear_form += equation


class FEMSolver:
    """Solves boundary value problem.

    Parameters
    ----------
    bilinear_form : ngsolve.BilinearForm
    linear_form : ngsolve.LinearForm
    solver : str
    preconditioner : str
    preconditioner_options: dict

    Notes
    -----

    The bilinear and linear form are
    used to assemble the linear system of equations to
    be solved by a (preconditioned) solver.
    The preconditioner options have to be chosen according
    to the chosen preconditioner and for now apply only for the
    `bddc` preconditioner.

    """

    def __init__(self,
                 bilinear_form: ngsolve.BilinearForm,
                 linear_form: ngsolve.LinearForm,
                 solver: str,
                 preconditioner: str,
                 preconditioner_options: str = {}) -> None:

        self.__a = bilinear_form
        self.__f = linear_form
        # TODO Physics related checks and asserts
        # For example: for floating not every preconditioner is suitable
        self.__preconditioner = ngsolve.Preconditioner(bf=self.__a,
                                                       type=preconditioner,
                                                       **preconditioner_options)
        if solver in ["CG", "GMRES"]:
            self.__solver = solver
        else:
            raise RuntimeError("The solver has to be either CG or GMRES")

    def solve_bvp(self, input: ngsolve.comp.GridFunction) \
            -> ngsolve.la.DynamicVectorExpression:
        """Solve the boundary value problem.

        Parameters
        ----------
        input : ngsolve.comp.GridFunction
            Solution vector with boundary values.

        Returns
        -------
        ngsolve.la.DynamicVectorExpression
            Numerical solution of the boundary value problem.
        """

        self.__a.Assemble()
        self.__f.Assemble()
        if self.__solver == "CG":
            # TODO make printrates, maxsteps and precision flexible => input dict
            inverse = ngsolve.CGSolver(mat=self.__a.mat,
                                       pre=self.__preconditioner.mat,
                                       printrates=True,
                                       maxsteps=10000,
                                       precision=1e-12)
        elif self.__solver == "GMRES":
            # TODO make printrates, maxsteps and precision flexible => input dict
            inverse = ngsolve.GMRESSolver(mat=self.__a.mat,
                                          pre=self.__preconditioner.mat,
                                          printrates=True,
                                          maxsteps=10000,
                                          precision=1e-12)
        r = self.__f.vec.CreateVector()
        r.data = self.__f.vec - self.__a.mat * input.vec
        return input.vec.data + inverse * r
