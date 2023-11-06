import logging
import os
import time
from abc import ABC, abstractmethod
from typing import List, Optional

import ngsolve
import numpy as np
import pandas as pd

from ossdbs.fem.mesh import Mesh
from ossdbs.fem.solver import Solver
from ossdbs.model_geometry import Contacts, ModelGeometry
from ossdbs.point_analysis import Lattice, Pathway, PointModel, VoxelLattice
from ossdbs.stimulation_signals import FrequencyDomainSignal
from ossdbs.utils.vtk_export import FieldSolution

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductor(ABC):
    """Template class of a volume conductor.

    Parameters
    ----------
    mesh : Mesh
    conductivity : ConductivityCF
    model_geometry : ModelGeometry
    solver : Solver

    # TODO add more abstractmethod ?
    """

    def __init__(
        self,
        geometry: ModelGeometry,
        conductivity: ConductivityCF,
        solver: Solver,
        order: int,
        meshing_parameters: dict,
        frequency_domain_signal: FrequencyDomainSignal,
    ) -> None:
        self._solver = solver
        self._model_geometry = geometry
        # to update contacts later, shall not be changed
        self._base_contacts = Contacts(geometry.contacts)
        # can be changed at every simulation run
        self._contacts = Contacts(geometry.contacts)
        self._order = order

        self._conductivity_cf = conductivity
        self._complex = conductivity.is_complex
        self.signal = frequency_domain_signal

        self._impedances = None
        self._mesh = Mesh(self._model_geometry.geometry, self._order)
        if meshing_parameters["LoadMesh"]:
            # TODO check that LoadPath is there if LoadMesh
            self._mesh.load_mesh(meshing_parameters["LoadPath"])
        else:
            self._mesh.generate_mesh(meshing_parameters["MeshingHypothesis"])

        if meshing_parameters["SaveMesh"]:
            self._mesh.save(meshing_parameters["SavePath"])

        # to save previous solution and do post-processing
        self._frequency = None
        self._sigma = None

        # to write results
        self._output_path = None

    @abstractmethod
    def compute_solution(self, frequency: float) -> None:
        """Compute solution at frequency."""
        pass

    def run_full_analysis(
        self,
        compute_impedance: bool = False,
        export_vtk: bool = False,
        point_model: PointModel = None,
        template_space: bool = False,
        activation_threshold: Optional[float] = None,
    ) -> dict:
        """Run entire volume conductor model.

        Notes
        -----
        TODO full documentation
        """
        timings: dict = {}
        if point_model is not None:
            timings["FieldExport"] = []
        if export_vtk:
            timings["VTKExport"] = []
        timings["ComputeSolution"] = []

        if compute_impedance:
            if self.is_complex:
                self._impedances = np.ndarray(
                    shape=(len(self.signal.frequencies)), dtype=complex
                )
            else:
                self._impedances = np.ndarray(shape=(len(self.signal.frequencies)))
        for idx, frequency in enumerate(self.signal.frequencies):
            _logger.info(f"Computing at frequency: {frequency}")
            if not self.current_controlled:
                _logger.debug("Get scaled voltage values")
                voltage_values = self.get_scaled_active_contact_voltages(
                    self.signal.amplitudes[idx]
                )
                self.update_contacts(voltages=voltage_values)
            else:
                if len(self.contacts.active) == 2:
                    for contact_idx, contact in enumerate(self.contacts.active):
                        self.contacts[contact.name].voltage = float(contact_idx)
                else:
                    _logger.debug("Get scaled current values")
                    for contact in self.contacts.active:
                        if not np.isclose(contact.voltage, 0):
                            raise ValueError(
                                "In multicontact current-controlled mode, only ground voltage (0V) can be set on active contacts!"
                            )
                    current_values = self.get_scaled_contact_currents(
                        self.signal.amplitudes[idx]
                    )
                    self.update_contacts(currents=current_values)
            time_0 = time.time()
            self.compute_solution(frequency)
            time_1 = time.time()
            timings["ComputeSolution"].append(time_1 - time_0)
            time_0 = time_1
            if compute_impedance:
                self._impedances[idx] = self.compute_impedance()
            # multicontact is covered by floating approach
            if self.current_controlled and len(self.contacts.active) == 2:
                impedance = (
                    self.impedances[idx]
                    if compute_impedance
                    else self.compute_impedance()
                )
                # use Ohm's law U = Z * I
                # and that the Fourier coefficient for the current is known
                amplitude = self.contacts.active[0].current
                # use positive current by construction
                if self.is_complex:
                    sign = np.sign(amplitude.real)
                else:
                    sign = np.sign(amplitude)
                amplitude *= sign

                scale_voltage = impedance * amplitude
                # directly access GridFunction here
                self._potential.vec[:] = (
                    scale_voltage * self._potential.vec.FV().NumPy()
                )
            if _logger.getEffectiveLevel() == logging.DEBUG:
                estimated_currents = self.estimate_currents()
                _logger.debug(
                    f"Estimated currents through contacts: {estimated_currents}"
                )
            if export_vtk:
                self.vtk_export()
                time_1 = time.time()
                timings["VTKExport"].append(time_1 - time_0)
                time_0 = time_1
            if point_model is not None:
                self.export_points(point_model, activation_threshold, template_space)
                time_1 = time.time()
                timings["FieldExport"] = time_1 - time_0
                time_0 = time_1

        # save impedance at all frequencies to file!
        if compute_impedance:
            _logger.info("Saving impedance")
            df = pd.DataFrame(
                {
                    "freq": self.signal.frequencies,
                    "real": self.impedances.real,
                    "imag": self.impedances.imag,
                }
            )
            df.to_csv(os.path.join(self.output_path, "impedance.csv"), index=False)

        return timings

    @property
    def output_path(self) -> str:
        """Returns the path to output."""
        return self._output_path

    @output_path.setter
    def output_path(self, path: str) -> None:
        """Set the path to write output.

        Notes
        -----
        Creates directory if it doesn't exist.
        """
        if not os.path.exists(path):
            os.makedirs(path)
        self._output_path = path

    @property
    def conductivity_cf(self) -> ConductivityCF:
        """Retruns the coefficient function of the conductivity."""
        return self._conductivity_cf

    @property
    def signal(self) -> FrequencyDomainSignal:
        return self._signal

    @signal.setter
    def signal(self, new_signal: FrequencyDomainSignal) -> None:
        self._check_signal(new_signal)
        self._signal = new_signal

    def _check_signal(self, new_signal: FrequencyDomainSignal) -> None:
        if new_signal.current_controlled:
            sum_currents = 0.0
            voltages_active = np.zeros(len(self.contacts.active))
            for idx, contact in enumerate(self.contacts.active):
                sum_currents += contact.current
                voltages_active[idx] = contact.voltage
            for contact in self.contacts.floating:
                if not np.all(np.isclose(voltages_active, 0.0)):
                    raise ValueError(
                        "In multipolar current-controlled mode, all active contacts have to be grounded!"
                    )
                sum_currents += contact.current
            if not np.isclose(sum_currents, 0):
                raise ValueError("The sum of all currents is not zero!")

    @property
    def current_controlled(self) -> bool:
        return self.signal.current_controlled

    @property
    def impedances(self) -> np.ndarray:
        return self._impedances

    @property
    def conductivity(self) -> ngsolve.CoefficientFunction:
        """Return conductivity of latest solution."""
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
            contact_voltages[contact.name] = (
                self._base_contacts[contact.name].voltage * factor
            )
        return contact_voltages

    def get_scaled_contact_currents(self, factor) -> dict:
        contact_currents = {}
        contacts = list(self.contacts.active)
        contacts_floating = list(self.contacts.floating)
        contacts.extend(contacts_floating)

        total_current = 0
        for contact in contacts:
            contact_current = self._base_contacts[contact.name].current * factor
            contact_currents[contact.name] = contact_current
            total_current += contact_current
        if not np.isclose(total_current, 0):
            raise RuntimeError("Sum of currents is not zero!")
        return contact_currents

    def update_contacts(
        self,
        voltages: Optional[dict] = None,
        currents: Optional[dict] = None,
        surface_impedances: Optional[dict] = None,
    ) -> None:
        """TODO document."""
        if surface_impedances is None:
            surface_impedances = {}
        if currents is None:
            currents = {}
        if voltages is None:
            voltages = {}
        if self.is_complex:
            self._contacts.voltages = voltages
            self._contacts.currents = currents
            self._contacts.surface_impedances = surface_impedances
            return
        self._contacts.voltages = np.real(voltages)
        self._contacts.currents = np.real(currents)
        self._contacts.surface_impedances = np.real(surface_impedances)

    @property
    def potential(self) -> ngsolve.GridFunction:
        return self._potential

    @property
    def frequency(self) -> float:
        return self._frequency

    def evaluate_potential_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric potential at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : Nx3 numpy.ndarray of lattice points

        Returns
        -------
        Nx1 numpy.ndarray

        Notes
        -----
        Make sure that points outside of the computational domain
        are filtered!
        """
        mesh = self.mesh.ngsolvemesh
        x, y, z = lattice.T
        pots = self.potential(mesh(x, y, z))
        return pots

    def evaluate_field_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric field components at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : Nx3 numpy.ndarray of lattice points

        Returns
        -------
        Nx3 numpy.ndarray

        TODO: filter points outside of the computational domain
        """
        mesh = self.mesh.ngsolvemesh
        x, y, z = lattice.T
        fields = self.electric_field(mesh(x, y, z))
        return fields

    @property
    def current_density(self) -> ngsolve.GridFunction:
        """Return current density in A/mm^2."""
        # scale to account for mm as length unit (not yet contained in conductivity)
        return 1e-3 * self.conductivity * self.electric_field

    @property
    def electric_field(self) -> ngsolve.GridFunction:
        return -ngsolve.grad(self.potential)

    def compute_impedance(self) -> complex:
        """Compute impedance at most recent solution.

        Notes
        -----
        The impedance is so far only available for two
        active contacts. It is computed by volume integration.
        This approach is superior to integration of the
        normal current density. It has been described for
        example in [Zimmermann2021a]_.
        Since the voltage drop is not known, we infer it
        from the voltages of the two contacts.
        By construction, the voltage is a positive value
        (or in the complex case, the real part).

        References
        ----------
        .. [Zimmermann2021a] Zimmermann, J., et al. (2021).
                             Frontiers in Bioengineering and Biotechnology, 9, 765516.
                             https://doi.org/10.3389/fbioe.2021.765516

        """
        if len(self.contacts.active) == 2:
            mesh = self._mesh.ngsolvemesh
            # do not need to account for mm because of integration
            power = ngsolve.Integrate(self.electric_field * self.current_density, mesh)
            # TODO integrate surface impedance by thin layer
            voltage = 0
            for idx, contact in enumerate(self.contacts.active):
                voltage += (-1) ** idx * contact.voltage
            if self.is_complex:
                sign = np.sign(voltage.real)
            else:
                sign = np.sign(voltage)
            return sign * voltage / power
        else:
            # TODO implement meaningful way to access contribution of individual
            # electrode to impedance
            raise NotImplementedError(
                "Impedance for more than two active contacts not yet supported"
            )

    def estimate_currents(self) -> dict:
        """Estimate currents by integration of normal component.

        Notes
        -----
        Meant for debugging purposes.
        If singularities are present, this method will not be accurate.
        """
        normal_vector = ngsolve.specialcf.normal(3)
        estimated_currents = {}
        for contact in self.contacts:
            # use that normal_vector always points outwards
            normal_current_density = -normal_vector * ngsolve.BoundaryFromVolumeCF(
                self.current_density
            )
            current = ngsolve.Integrate(
                normal_current_density * ngsolve.ds(contact.name), self.mesh.ngsolvemesh
            )
            estimated_currents[contact.name] = current
        return estimated_currents

    def vtk_export(self) -> None:
        """Export all relevant properties to VTK."""
        ngmesh = self.mesh.ngsolvemesh
        # TODO add frequency to name
        FieldSolution(self.potential, "potential", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "potential")
        )

        FieldSolution(self.electric_field, "E-field", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "E-field")
        )

        if self.conductivity_cf.is_tensor:
            # Naming convention by ParaView!
            cf_list = (
                self.conductivity[0],  # xx
                self.conductivity[4],  # yy
                self.conductivity[8],  # zz
                self.conductivity[1],  # xy
                self.conductivity[5],  # yz
                self.conductivity[2],  # xz
            )
            conductivity_export = ngsolve.CoefficientFunction(cf_list, dims=(6,))
        else:
            conductivity_export = self.conductivity
        FieldSolution(
            conductivity_export, "conductivity", ngmesh, self.is_complex
        ).save(os.path.join(self.output_path, "conductivity"))

        if self.conductivity_cf.is_tensor:
            dti_voxel = self.conductivity_cf.dti_voxel_distribution
            # Naming convention by ParaView!
            cf_list = (
                dti_voxel[0],  # xx
                dti_voxel[4],  # yy
                dti_voxel[8],  # zz
                dti_voxel[1],  # xy
                dti_voxel[5],  # yz
                dti_voxel[2],  # xz
            )
            dti_export = ngsolve.CoefficientFunction(cf_list, dims=(6,))
            FieldSolution(dti_export, "dti", ngmesh, False).save(
                os.path.join(self.output_path, "dti")
            )

        FieldSolution(
            self.conductivity_cf.material_distribution(self.mesh),
            "material",
            ngmesh,
            False,
        ).save(os.path.join(self.output_path, "material"))

    def export_points(
        self, point_model, activation_threshold: float, template_space: bool
    ) -> None:
        """Export all relevant properties to CSV, HDF5, or Nifty format.

        Notes
        -----
        Only for pathways a HDF5 file is generated and only for the lattice
        model, a Nifty file is generated

        Parameters
        ----------
        point_model: oss_dbs point_model
            Contains the points on which to evaluate the solution

        activation_threshold: float

        template_space: bool
        """
        grid_pts = point_model.points_in_mesh(self.mesh)
        lattice_mask = np.invert(grid_pts.mask)
        lattice = point_model.filter_for_geometry(grid_pts)
        inside_csf = self.get_points_in_csf(lattice)
        inside_encap = self.get_points_in_encapsulation_layer(lattice)
        if isinstance(point_model, Pathway):
            # if pathway, always mark complete axons
            axon_mask = point_model._axon_mask
            [inside_csf, inside_encap] = point_model.filter_csf_encap(
                inside_csf, inside_encap
            )
            # create index for axons
            index = point_model.create_index(lattice)
        else:
            index = np.reshape(np.arange(len(lattice)), (len(lattice), 1))

        potentials = self.evaluate_potential_at_points(lattice)

        fields = self.evaluate_field_at_points(lattice)
        field_mags = np.linalg.norm(fields, axis=1).reshape((fields.shape[0], 1))

        # CSV exports
        df_pot = pd.DataFrame(
            np.concatenate(
                [
                    index,
                    lattice,
                    potentials.reshape((potentials.shape[0], 1)),
                    inside_csf,
                    inside_encap,
                ],
                axis=1,
            ),
            columns=[
                "index",
                "x-pt",
                "y-pt",
                "z-pt",
                "potential",
                "inside_csf",
                "inside_encap",
            ],
        )
        df_pot.to_csv(os.path.join(self.output_path, "oss_potentials.csv"), index=False)
        df_field = pd.DataFrame(
            np.concatenate(
                [index, lattice, fields, field_mags, inside_csf, inside_encap],
                axis=1,
            ),
            columns=[
                "index",
                "x-pt",
                "y-pt",
                "z-pt",
                "x-field",
                "y-field",
                "z-field",
                "magnitude",
                "inside_csf",
                "inside_encap",
            ],
        )
        if template_space:
            df_field.to_csv(
                os.path.join(self.output_path, "E_field_Template_space.csv"),
                index=False,
            )
        else:
            df_field.to_csv(
                os.path.join(self.output_path, "E_field_MRI_space.csv"),
                index=False,
            )

        # HDF5 exports only for Pathways
        if isinstance(point_model, Pathway):
            if axon_mask is None:
                raise RuntimeError("The creation of the axon_mask did not work")
            else:
                point_model.save_hdf5(
                    axon_mask,
                    lattice,
                    potentials,
                    fields,
                    field_mags,
                    self.output_path,
                )

        # Nifti exports
        if isinstance(point_model, VoxelLattice) or isinstance(point_model, Lattice):
            field_mags_full = np.zeros(lattice_mask.shape[0], float)
            # TODO set all values to -1?
            # field_mags_full -= 1
            # overwrite values inside mesh
            field_mags_full[lattice_mask[:, 0]] = field_mags[:, 0]

            if type(point_model) == VoxelLattice:
                suffix = ""
            elif isinstance(point_model, Lattice):
                suffix = "_WA"
            point_model.save_as_nifti(
                field_mags_full,
                os.path.join(self.output_path, f"E_field_solution{suffix}.nii"),
            )
            point_model.save_as_nifti(
                field_mags_full,
                os.path.join(self.output_path, f"VTA_solution{suffix}.nii"),
                binarize=True,
                activation_threshold=activation_threshold,
            )

    def floating_values(self) -> dict:
        """Read out floating potentials."""
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
        return ngsolve.HDiv(
            mesh=self.mesh, order=max(1, self._order - 1), complex=self._complex
        )

    def h1_space(self, boundaries: List[str]) -> ngsolve.H1:
        """Return a h1 space on the mesh.

        Parameters
        ----------
        boundaries : list of str
            List of boundary names.

        Returns
        -------
        ngsolve.H1
        """
        dirichlet = "|".join(boundary for boundary in boundaries)
        return ngsolve.H1(
            mesh=self.mesh.ngsolvemesh,
            order=self._order,
            dirichlet=dirichlet,
            complex=self._complex,
            wb_withedges=False,
        )

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space on the mesh.

        Returns
        -------
        ngsolve.NumberSpace
            Space with only one single (global) DOF.
        """
        return ngsolve.NumberSpace(
            mesh=self.mesh.ngsolvemesh, order=0, complex=self.is_complex
        )

    def get_points_in_encapsulation_layer(self, points: np.ndarray) -> np.ndarray:
        """Return mask for points in encapsulation layer."""
        encap_cf = self.mesh.ngsolvemesh.RegionCF(
            ngsolve.VOL, {"EncapsulationLayer_*": 1.0}, default=0
        )
        mesh = self.mesh.ngsolvemesh
        x, y, z = points.T
        return np.isclose(encap_cf(mesh(x, y, z)), 1.0)

    def get_points_in_csf(self, points: np.ndarray) -> np.ndarray:
        """Return mask for points in CSF."""
        material_distribution = self.conductivity_cf.material_distribution(self.mesh)
        mesh = self.mesh.ngsolvemesh
        x, y, z = points.T
        return np.isclose(
            material_distribution(mesh(x, y, z)), self.conductivity_cf.materials["CSF"]
        )
