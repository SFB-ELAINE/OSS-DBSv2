from abc import ABC, abstractmethod
import ngsolve
from ossdbs.model_geometry import ModelGeometry, Contacts
from ossdbs.fem.solver import Solver
from ossdbs.fem.mesh import Mesh
from ossdbs.stimulation_signals import FrequencyDomainSignal
from .conductivity import ConductivityCF
from typing import List
import numpy as np

import logging
_logger = logging.getLogger(__name__)


class VolumeConductor(ABC):
    """Template class of a volume conductor

    Parameters
    ----------
    mesh : Mesh
    conductivity : ConductivityCF
    model_geometry : ModelGeometry
    solver : Solver

    # TODO add more abstractmethod ?
    """

    def __init__(self,
                 geometry: ModelGeometry,
                 conductivity: ConductivityCF,
                 solver: Solver,
                 order: int,
                 meshing_parameters: dict,
                 frequency_domain_signal: FrequencyDomainSignal) -> None:
        self._solver = solver
        self._model_geometry = geometry
        # to update contacts later, shall not be changed
        self._base_contacts = Contacts(geometry.contacts)
        # can be changed at every simulation run
        self._contacts = Contacts(geometry.contacts)
        self._order = order

        self._conductivity_cf = conductivity
        self._complex = conductivity.is_complex
        self._signal = frequency_domain_signal

        self._impedances = None
        self._mesh = Mesh(self._model_geometry.geometry, self._order)
        if meshing_parameters["LoadMesh"]:
            # TODO check that LoadPath is there if LoadMesh
            self._mesh.load_mesh(meshing_parameters["LoadPath"])
        else:
            self._mesh.generate_mesh(meshing_parameters['MeshingHypothesis'])

        if meshing_parameters["SaveMesh"]:
            self._mesh.save_mesh(Mesh["Path"])

        # to save previous solution and do post-processing
        self._frequency = None
        self._sigma = None

    @abstractmethod
    def compute_solution(self,
                         frequency: float
                         ) -> None:
        pass

    def run_full_analysis(self,
                          compute_impedance: bool = False,
                          lattice: np.ndarray = None
                          ) -> None:
        if compute_impedance:
            self._impedances = np.ndarray(shape=(len(self.signal.frequencies),))
        for idx, frequency in enumerate(self.signal.frequencies):
            _logger.info("Computing at frequency: {}".format(frequency))
            if not self.current_controlled:
                voltage_values = self.get_scaled_active_contact_voltages(self.signal.amplitudes[idx])
                self.update_contacts(voltages=voltage_values)
            self.compute_solution(frequency)

            if self.current_controlled:
                # TODO implement current control
                # Steps: compute impedance and rescale
                raise NotImplementedError("Current-controlled mode not yet implemented")
                continue
            if compute_impedance:
                self._impedances[idx] = self.compute_impedance()

    @property
    def conductivity_cf(self) -> ConductivityCF:
        return self._conductivity_cf

    @property
    def signal(self) -> FrequencyDomainSignal:
        return self._signal

    @property
    def current_controlled(self) -> bool:
        return self.signal.current_controlled

    @property
    def impedances(self) -> np.ndarray:
        return self._impedances

    @property
    def conductivity(self) -> ngsolve.CoefficientFunction:
        """Return conductivity of latest solution
        """
        return self._sigma

    @property
    def is_complex(self) -> bool:
        """Return the state of the data type for spaces. True if complex,
        False otherwise.

        Returns
        -------
        bool
        """
        return self._complex

    @is_complex.setter
    def is_complex(self, value: bool) -> None:
        self._complex = value

    @property
    def model_geometry(self) -> ModelGeometry:
        return self._model_geometry

    @property
    def mesh(self) -> Mesh:
        return self._mesh

    @property
    def solver(self) -> Solver:
        return self._solver

    @property
    def contacts(self) -> Contacts:
        return self._contacts

    def get_scaled_active_contact_voltages(self, factor) -> dict:
        contact_voltages = {}
        for contact in self.contacts.active:
            contact_voltages[contact.name] = self._base_contacts[contact.name].voltage * factor
        return contact_voltages

    def update_contacts(self, voltages: dict = {}, currents: dict = {}, surface_impedances: dict = {}) -> None:
        """TODO document
        """
        self._contacts.voltages = voltages
        self._contacts.currents = currents
        self._contacts.surface_impedances = surface_impedances

    @property
    def potential(self) -> ngsolve.GridFunction:
        return self._potential

    @property
    def frequency(self) -> float:
        return self._frequency

    def evaluate_potential_at_points(self, lattice: np.ndarray) -> np.ndarray:
        # TODO
        pass

    @property
    def current_density(self) -> ngsolve.GridFunction:
        _logger.debug("Compute current density at frequency {} Hz".format(self.frequency))
        return self.conductivity * self.electric_field

    @property
    def electric_field(self) -> ngsolve.GridFunction:
        return -ngsolve.grad(self.potential)

    def compute_impedance(self) -> complex:
        """TODO document
        """
        if len(self.contacts.active) == 2:
            curr_dens_conj = ngsolve.Conj(self.current_density)
            mesh = self._mesh.ngsolvemesh
            power = ngsolve.Integrate(self.electric_field * curr_dens_conj, mesh)
            # TODO integrate surface impedance
            voltage = 0
            for idx, contact in enumerate(self.contacts.active):
                voltage += (-1)**idx * contact.voltage
            return voltage / power
        else:
            raise NotImplementedError("Impedance for more than two active contacts not yet supported")

    def floating_values(self) -> dict:
        floating_voltages = {}
        for contact in self.contacts.floating:
            floating_voltages[contact.name] = contact.voltage
        return floating_voltages

    def flux_space(self) -> ngsolve.comp.HDiv:
        """Return a flux space on the mesh.

        Returns
        -------
        ngsolve.HDiv

        Notes
        -----

        The HDiv space is returned with a minimum order of 1.
        It is needed for the a-posteriori error estimator
        needed for adaptive mesh refinement.

        """

        return ngsolve.HDiv(mesh=self.mesh,
                            order=max(1, self._order - 1),
                            complex=self._complex)

    def h1_space(self, boundaries: List[str]) -> ngsolve.comp.H1:
        """Return a h1 space on the mesh.

        Parameters
        ----------
        boundaries : list of str
            List of boundary names.

        Returns
        -------
        ngsolve.H1
        """

        dirichlet = '|'.join(boundary for boundary in boundaries)
        return ngsolve.H1(mesh=self.mesh.ngsolvemesh,
                          order=self._order,
                          dirichlet=dirichlet,
                          complex=self._complex,
                          wb_withedges=False)

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space on the mesh.

        Returns
        -------
        ngsolve.NumberSpace
            Space with only one single (global) DOF.
        """

        return ngsolve.NumberSpace(mesh=self.mesh.ngsolvemesh,
                                   order=0,
                                   complex=self.is_complex)

    def surfacel2_space(self, boundaries: List[str]) -> ngsolve.comp.SurfaceL2:
        """Return a number SurfaceL2 on the mesh.

        Returns
        -------
        ngsolve.SurfaceL2

        Notes
        -----

        The SurfaceL2 space is returned with a minimum order of 1.
        It is needed to impose floating potentials.

        """

        dirichlet = '|'.join(boundary for boundary in boundaries)
        return ngsolve.SurfaceL2(mesh=self.mesh.ngsolvemesh,
                                 order=max(1, self._order - 1),
                                 dirichlet=dirichlet,
                                 complex=self.is_complex)
