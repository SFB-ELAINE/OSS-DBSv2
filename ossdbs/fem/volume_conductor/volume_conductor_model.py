import logging
import os
import time
from abc import ABC, abstractmethod
from typing import List, Optional, Tuple, Union

import ngsolve
import numpy as np
import pandas as pd
import scipy.fft as fft

from ossdbs.fem.mesh import Mesh
from ossdbs.fem.solver import Solver
from ossdbs.model_geometry import Contacts, ModelGeometry
from ossdbs.point_analysis import Lattice, Pathway, PointModel, TimeResult, VoxelLattice
from ossdbs.stimulation_signals import (
    FrequencyDomainSignal,
    get_indices_in_octave_band,
    get_octave_band_indices,
    get_timesteps,
)
from ossdbs.utils.collapse_vta import get_collapsed_VTA
from ossdbs.utils.vtk_export import FieldSolution

from .conductivity import ConductivityCF

_logger = logging.getLogger(__name__)


class VolumeConductor(ABC):
    """Template class of a volume conductor.

    Parameters
    ----------
    geometry : ModelGeometry
        Model geometry, brain with implanted electrodes
    conductivity : ConductivityCF
        Material information
    solver : Solver
        Solver (linear algebra part)
    order: int
        Order of solver and mesh (curved elements)
    meshing_parameters: dict
        Dictionary with setting for meshing
    frequency_domain_signal: FrequencyDomainSignal
        Representation of signal in frequency domain
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
        self._order = order

        self._model_geometry = geometry

        self._signal = frequency_domain_signal
        _logger.debug(f"Signal with amplitudes: {self.signal.amplitudes}")
        # contacts of electrode
        self._contacts = Contacts(geometry.contacts)
        # set voltages for two-contact current controlled stimulation
        # check if current-controlled mode is correctly prepared
        if self.current_controlled:
            self.prepare_current_controlled_mode()

        _logger.debug(f"Assigned base contacts with properties:\n {self._contacts}")

        self._conductivity_cf = conductivity
        self._complex = conductivity.is_complex

        # to store impedances at all frequencies
        self._impedances = None
        # to store the voltage (current-controlled)
        # or current (voltage-controlled) stimulation
        self._free_stimulation_variable = None
        self._stimulation_variable = None

        # generate the mesh
        self._mesh = Mesh(self._model_geometry.geometry, self._order)
        if meshing_parameters["LoadMesh"]:
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

        # frequency at which VTK shall be exported
        self._export_frequency = None

    @abstractmethod
    def compute_solution(self, frequency: float) -> None:
        """Compute solution at frequency.

        Parameters
        ----------
        frequency: float
            Frequency at which solution is computed
        """
        pass

    def run_full_analysis(
        self,
        compute_impedance: bool = False,
        export_vtk: bool = False,
        point_model: PointModel = None,
        activation_threshold: Optional[float] = None,
        collapse_vta: bool = False,
    ) -> dict:
        """Run volume conductor model at all frequencies.

        Parameters
        ----------
        compute_impedance: bool
            If True, the impedance will be computed at each frequency.
        export_vtk: bool
            VTK export for visualization in ParaView
        point_model: PointModel
            PointModel to extract solution at points for VTA / PAM
        activation_threshold: float
            If VTA is estimated by threshold, provide it here.
        collapse_vta: bool
            Whether to collapse VTA

        Notes
        -----
        TODO full documentation
        TODO pass list of PointModel to have more than one evaluated (PAM + VTA)
        The volume conductor model is run at all frequencies
        and the time-domain signal is computed (if relevant).
        """
        timings = self.setup_timings_dict(export_vtk, point_model)

        dtype = float
        if self.is_complex:
            dtype = complex

        # always compute impedance for CC with 2 contacts
        if self.current_controlled and len(self.contacts.active) == 2:
            _logger.info(
                """Set compute_impedance to True.
                   Impedance calculation is required for 2 contacts."""
            )
            compute_impedance = True

        if self.signal.octave_band_approximation:
            frequency_indices = get_octave_band_indices(self.signal.frequencies)
            # add DC component
            if not np.isclose(self.signal.amplitudes[0], 0.0):
                frequency_indices = np.insert(frequency_indices, 0, 0)
        else:
            frequency_indices = np.arange(len(self.signal.frequencies))

        median_index = np.median(frequency_indices)
        self._export_frequency = self.signal.frequencies[int(median_index)]

        self._free_stimulation_variable = np.zeros(
            shape=(len(self.signal.frequencies), len(self.contacts.active)),
            dtype=complex,
        )
        self._stimulation_variable = np.zeros(
            shape=(len(self.signal.frequencies), len(self.contacts.active)),
            dtype=complex,
        )

        if compute_impedance:
            self._impedances = np.ndarray(
                shape=(len(self.signal.frequencies)), dtype=dtype
            )

        # TODO wrap in function
        if point_model is not None:
            (
                grid_pts,
                lattice_mask,
                lattice,
                inside_csf,
                inside_encap,
                axon_index,
            ) = self.prepare_point_model_grids(point_model)
            export_collapsed_vta = isinstance(point_model, Lattice) and collapse_vta
            export_nifti = isinstance(point_model, VoxelLattice) or isinstance(
                point_model, Lattice
            )
            tmp_potential_freq_domain = np.zeros(
                shape=(len(self.signal.frequencies), len(lattice)), dtype=complex
            )
            tmp_Ex_freq_domain = np.zeros(
                shape=(len(self.signal.frequencies), len(lattice)), dtype=complex
            )
            tmp_Ey_freq_domain = np.zeros(
                shape=(len(self.signal.frequencies), len(lattice)), dtype=complex
            )
            tmp_Ez_freq_domain = np.zeros(
                shape=(len(self.signal.frequencies), len(lattice)), dtype=complex
            )

        for freq_idx in frequency_indices:
            frequency = self.signal.frequencies[freq_idx]
            _logger.info(f"Computing at frequency: {frequency}")
            time_0 = time.time()

            # prepare storing at multiple frequencies
            if self.signal.octave_band_approximation:
                band_indices = get_indices_in_octave_band(
                    freq_idx, frequency_indices, len(self.signal.frequencies) - 1
                )
                _logger.debug(
                    f"""Band indices from {band_indices[0]}
                        to {band_indices[-1]}"""
                )
                _logger.debug(
                    f"""Band frequencies from {self.signal.frequencies[band_indices[0]]}
                        to {self.signal.frequencies[band_indices[-1]]}"""
                )

            else:
                band_indices = [freq_idx]

            # check if conductivity has changed
            sigma_has_changed = self._has_sigma_changed(freq_idx)
            if sigma_has_changed:
                self.compute_solution(frequency)
                if compute_impedance:
                    impedance = self.compute_impedance()
                    self._impedances[band_indices] = impedance
            else:
                _logger.info(f"Skipped computation at {frequency} Hz")
                if compute_impedance:
                    impedance = self._impedances[freq_idx - 1]
                    self._impedances[band_indices] = impedance

            # scale factor: is one for VC and depends on impedance for other case
            self._scale_factor = self.get_scale_factor(freq_idx)
            _logger.debug(f"Scale factor: {self._scale_factor}")

            self._store_solution_at_contacts(band_indices)
            if _logger.getEffectiveLevel() == logging.DEBUG:
                estimated_currents = self.estimate_currents()
                _logger.debug(
                    f"Estimated currents through contacts: {estimated_currents}"
                )

            time_1 = time.time()
            timings["ComputeSolution"].append(time_1 - time_0)
            time_0 = time_1

            # exports in frequency domain
            if export_vtk and np.isclose(frequency, self._export_frequency):
                # TODO not yet scaled
                self.vtk_export()
                time_1 = time.time()
                timings["VTKExport"].append(time_1 - time_0)
                time_0 = time_1

            # export point results
            if point_model is not None:
                potentials = self.evaluate_potential_at_points(lattice)
                fields = self.evaluate_field_at_points(lattice)
                field_mags = self.compute_field_magnitude(fields)

                # copy values for time-domain analysis
                self._copy_values_for_time_domain(
                    band_indices, tmp_potential_freq_domain, potentials[:, 0]
                )
                self._copy_values_for_time_domain(
                    band_indices, tmp_Ex_freq_domain, fields[:, 0]
                )
                self._copy_values_for_time_domain(
                    band_indices, tmp_Ey_freq_domain, fields[:, 1]
                )
                self._copy_values_for_time_domain(
                    band_indices, tmp_Ez_freq_domain, fields[:, 2]
                )

                time_1 = time.time()
                timings["CopyValues"].append(time_1 - time_0)
                time_0 = time_1

                # export results: VTA
                if np.isclose(frequency, self._export_frequency):
                    self.export_potential_to_csv(
                        frequency,
                        potentials,
                        axon_index,
                        lattice,
                        inside_csf,
                        inside_encap,
                    )
                    self.export_field_to_csv(
                        frequency,
                        fields,
                        field_mags,
                        axon_index,
                        lattice,
                        inside_csf,
                        inside_encap,
                        export_collapsed_vta,
                    )
                    if export_nifti:
                        self.export_nifti_files(
                            field_mags, point_model, activation_threshold, lattice_mask
                        )
                    # HDF5 exports only for Pathways
                    if isinstance(point_model, Pathway):
                        self.export_to_hdf5(
                            point_model, lattice, potentials, fields, field_mags
                        )
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

        if point_model and len(self.signal.frequencies) > 1:
            _logger.info("Reconstructing time-domain signal.")
            timesteps = get_timesteps(
                self.signal.cutoff_frequency,
                self.signal.base_frequency,
                self.signal.signal_length,
            )
            # TODO check if re-ordering is possible
            potential_in_time = self.reconstruct_time_signals(
                len(lattice), tmp_potential_freq_domain
            )
            potential_in_time = np.swapaxes(potential_in_time, 0, 1)
            Ex_in_time = self.reconstruct_time_signals(len(lattice), tmp_Ex_freq_domain)
            Ex_in_time = np.swapaxes(Ex_in_time, 0, 1)
            Ey_in_time = self.reconstruct_time_signals(len(lattice), tmp_Ey_freq_domain)
            Ey_in_time = np.swapaxes(Ey_in_time, 0, 1)
            Ez_in_time = self.reconstruct_time_signals(len(lattice), tmp_Ez_freq_domain)
            Ez_in_time = np.swapaxes(Ez_in_time, 0, 1)
            field_in_time = np.column_stack((Ex_in_time, Ey_in_time, Ez_in_time))
            self.create_time_result(
                point_model,
                lattice,
                timesteps,
                potential_in_time,
                field_in_time,
                inside_csf,
                inside_encap,
            )
            _logger.debug(f"Shapes of timesteps: {timesteps.shape}")
            _logger.debug(f"Shapes of potential_in_time: {potential_in_time.shape}")
            _logger.debug(f"Shapes of Ex_in_time: {Ex_in_time.shape}")
            _logger.debug(f"Shapes of Ey_in_time: {Ey_in_time.shape}")
            _logger.debug(f"Shapes of Ez_in_time: {Ez_in_time.shape}")
            _logger.debug(f"Shapes of field_in_time: {field_in_time.shape}")
            time_1 = time.time()
            timings["ReconstructTimeSignals"] = time_1 - time_0
            time_0 = time_1

        if len(self.signal.frequencies) > 1:
            self.export_solution_at_contacts()
        return timings

    def export_solution_at_contacts(self) -> None:
        """Time-domain solution export."""
        timesteps = get_timesteps(
            self.signal.cutoff_frequency,
            self.signal.base_frequency,
            self.signal.signal_length,
        )
        free_stimulation_variable_in_time = self.reconstruct_time_signals(
            len(self.contacts.active), self._free_stimulation_variable
        )
        stimulation_variable_in_time = self.reconstruct_time_signals(
            len(self.contacts.active), self._stimulation_variable
        )
        free_stimulation_variable_at_contact = {}
        free_stimulation_variable_at_contact["time"] = timesteps
        for contact_idx, contact in enumerate(self.contacts.active):
            free_stimulation_variable_at_contact[
                contact.name + "_free"
            ] = free_stimulation_variable_in_time[:, contact_idx]
            free_stimulation_variable_at_contact[
                contact.name
            ] = stimulation_variable_in_time[:, contact_idx]
        df = pd.DataFrame(free_stimulation_variable_at_contact)
        df.to_csv(
            os.path.join(self.output_path, "stimulation_in_time.csv"), index=False
        )

    def compute_field_magnitude(self, fields):
        """Compute magnitude of field vector."""
        if len(fields.shape) != 2 and fields.shape[1] != 3:
            raise ValueError(
                "Attempt to compute field magnitude from invalid field vector"
            )
        # TODO fix for complex vectors
        field_mags = np.linalg.norm(fields, axis=1).reshape((fields.shape[0], 1))
        _logger.debug(f"Shapes of field_mags: {field_mags.shape}")
        return field_mags

    def create_time_result(
        self,
        point_model,
        lattice,
        timesteps,
        potential_in_time,
        field_in_time,
        inside_csf,
        inside_encap,
    ) -> TimeResult:
        """Prepare time rsult and save it to file."""
        time_result = TimeResult(
            time_steps=timesteps,
            points=lattice,
            potential=potential_in_time,
            electric_field_vector=field_in_time,
            electric_field_magnitude=self.compute_field_magnitude(field_in_time),
            inside_csf=inside_csf,
            inside_encap=inside_encap,
        )
        point_model.save(
            time_result, os.path.join(self.output_path, "oss_time_result.h5")
        )
        _logger.info("Created time result and saved to file")

    def reconstruct_time_signals(
        self,
        n_lattice_points: int,
        freq_domain_signal: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute time signals from frequency-domain data.

        Notes
        -----
        If octave band approximation is used,
        it will also be copied to frequencies at other frequencies in band.

        """
        frequency_indices = self.signal.frequencies / self.signal.base_frequency
        frequency_indices = frequency_indices.astype(np.uint16)

        if n_lattice_points == 0:
            raise ValueError("No lattice points to process.")

        # For even signals, highest frequency is not in positive frequencies
        positive_freqs_part = freq_domain_signal
        if self.signal.signal_length % 2 == 0:
            positive_freqs_part = freq_domain_signal[:-1]

        # Append the reverted signal without the DC frequency
        tmp_freq_domain = np.append(
            positive_freqs_part,
            np.conjugate(np.flip(freq_domain_signal[1:], axis=0)),
            axis=0,
        )
        # run ifft with maximum possible amount of workers
        result_in_time = fft.ifft(tmp_freq_domain, axis=0, workers=-1)
        return result_in_time.real

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
        """Returns the coefficient function of the conductivity."""
        return self._conductivity_cf

    @property
    def signal(self) -> FrequencyDomainSignal:
        """Returns the frequency-domain representation of stimulation signal."""
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
                all_active_contacts_grounded = np.all(np.isclose(voltages_active, 0.0))
                if not all_active_contacts_grounded:
                    raise ValueError(
                        """In multipolar current-controlled mode,
                           all active contacts have to be grounded!"""
                    )
                sum_currents += contact.current
            if not np.isclose(sum_currents, 0):
                raise ValueError("The sum of all currents is not zero!")

    @property
    def current_controlled(self) -> bool:
        """Return if stimulation is current-controlled."""
        return self.signal.current_controlled

    @property
    def impedances(self) -> np.ndarray:
        """Return list of impedances."""
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
        """If complex mode (EQS) is used or not."""
        self._complex = value

    @property
    def model_geometry(self) -> ModelGeometry:
        """The underlying model geometry used for mesh generation."""
        return self._model_geometry

    @property
    def mesh(self) -> Mesh:
        """The mesh used in computations."""
        return self._mesh

    @property
    def solver(self) -> Solver:
        """The solver used in the VCM."""
        return self._solver

    @property
    def contacts(self) -> Contacts:
        """A list of contacts in the VCM."""
        return self._contacts

    @property
    def potential(self) -> ngsolve.GridFunction:
        """Return solution at most recent frequency."""
        return self._potential

    @property
    def frequency(self) -> float:
        """Most recent frequency, not equal to the frequency of the signal!."""
        return self._frequency

    def evaluate_potential_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric potential at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : np.ndarray
            Nx3 numpy.ndarray of lattice points

        Notes
        -----
        Requires that points outside of the computational domain
        have been filtered!
        """
        mesh = self.mesh.ngsolvemesh
        x, y, z = lattice.T
        pots = self.potential(mesh(x, y, z))
        return pots

    def evaluate_field_at_points(self, lattice: np.ndarray) -> np.ndarray:
        """Return electric field components at specifed 3-D coordinates.

        Parameters
        ----------
        lattice : np.ndarray
            Nx3 numpy.ndarray of lattice points

        Notes
        -----
        Requires that points outside of the computational domain
        have been filtered!
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
        """Compute electric field from potential."""
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
        FieldSolution(self.potential, "potential", ngmesh, self.is_complex).save(
            os.path.join(self.output_path, "potential")
        )

        FieldSolution(self.electric_field, "E_field", ngmesh, self.is_complex).save(
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

    def prepare_point_model_grids(self, point_model):
        """Prepare grid for point evaluation and label points outside grid.

        Parameters
        ----------
        point_model: PointModel
            PointModel to be used

        Returns
        -------
        grid_pts: TODO documentation
        """
        grid_pts = point_model.points_in_mesh(self.mesh)
        lattice_mask = np.invert(grid_pts.mask)
        lattice = point_model.filter_for_geometry(grid_pts)
        inside_csf = self.get_points_in_csf(lattice)
        inside_encap = self.get_points_in_encapsulation_layer(lattice)
        if isinstance(point_model, Pathway):
            # if pathway, always mark complete axons
            point_model.filter_csf_encap(inside_csf, inside_encap)
            # create index for axons
            axon_index = point_model.create_index(lattice)
        else:
            axon_index = np.reshape(np.arange(len(lattice)), (len(lattice), 1))
        return grid_pts, lattice_mask, lattice, inside_csf, inside_encap, axon_index

    def export_potential_to_csv(
        self,
        frequency: float,
        potentials: np.ndarray,
        axon_index: np.ndarray,
        lattice: np.ndarray,
        inside_csf: np.ndarray,
        inside_encap: np.ndarray,
    ) -> None:
        """Export potential to CSV.

        Notes
        -----
        Only for pathways a HDF5 file is generated and only for the lattice
        model, a Nifti file is generated
        TODO Documentation needed!!!

        Parameters
        ----------
        frequency: float
            Frequency of exported solution
        potentials: np.ndarray
            List of potentials
        axon_index: np.ndarray
            TODO
        lattice: np.ndarray:
            TODO
        inside_csf: np.ndarray:
            TODO
        inside_encap: np.ndarray:
            TODO
        """
        df_pot = pd.DataFrame(
            np.concatenate(
                [
                    axon_index,
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
        # save frequency
        df_pot["frequency"] = frequency
        df_pot.to_csv(os.path.join(self.output_path, "oss_potentials.csv"), index=False)

    def export_field_to_csv(
        self,
        frequency: float,
        fields: np.ndarray,
        field_mags: np.ndarray,
        axon_index: np.ndarray,
        lattice: np.ndarray,
        inside_csf: np.ndarray,
        inside_encap: np.ndarray,
        collapse_VTA: bool,
    ):
        """Write field values to CSV.

        Parameters
        ----------
        frequency: float
            Frequency of exported solution
        fields: np.ndarray
            List of field vectors
        field_mags: np.ndarray
            List of field magnitudes
        axon_index: np.ndarray
            TODO
        lattice: np.ndarray:
            TODO
        inside_csf: np.ndarray:
            TODO
        inside_encap: np.ndarray:
            TODO
        collapse_VTA: bool
            Whether to collapse VTA
        """
        df_field = pd.DataFrame(
            np.concatenate(
                [axon_index, lattice, fields, field_mags, inside_csf, inside_encap],
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
        # save frequency
        df_field["frequency"] = frequency

        if collapse_VTA:
            _logger.info("Collapse VTA by virtually removing the electrode")
            field_on_probed_points = np.concatenate(
                [lattice, fields, field_mags], axis=1
            )

            electrode = self.model_geometry.electrodes[0]
            implantation_coordinate = electrode._position
            lead_direction = electrode._direction
            lead_diam = electrode._parameters.lead_diameter

            field_on_probed_points_collapsed = get_collapsed_VTA(
                field_on_probed_points,
                implantation_coordinate,
                lead_direction,
                lead_diam,
            )

            df_collapsed_field = pd.DataFrame(
                np.concatenate(
                    [
                        axon_index,
                        field_on_probed_points_collapsed,
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
                    "x-field",
                    "y-field",
                    "z-field",
                    "magnitude",
                    "inside_csf",
                    "inside_encap",
                ],
            )
            df_collapsed_field.to_csv(
                os.path.join(self.output_path, "E_field.csv"),
                index=False,
            )
        else:
            df_field.to_csv(
                os.path.join(self.output_path, "E_field.csv"),
                index=False,
            )

    def export_to_hdf5(self, point_model, lattice, potentials, fields, field_mags):
        """Export solution to HDF5 file."""
        point_model.save_hdf5(lattice, potentials, fields, field_mags, self.output_path)

    def export_nifti_files(
        self,
        field_mags: np.ndarray,
        point_model: PointModel,
        activation_threshold: float,
        lattice_mask: np.ndarray,
    ) -> None:
        """Export field to Nifti format.

        Parameters
        ----------
        field_mags: np.ndarray
            Field magnitudes to export
        point_model: oss_dbs point_model
            Contains the points on which to evaluate the solution
        activation_threshold: float
            Threshold to define VTA
        lattice_mask: np.ndarray
            Mask points outside domain

        activation_threshold: float
        """
        field_mags_full = np.zeros(lattice_mask.shape[0])
        field_mags_full[lattice_mask[:, 0]] = np.real(field_mags[:, 0]) * 1000.0

        point_model.save_as_nifti(
            field_mags_full,
            os.path.join(self.output_path, "E_field_solution.nii"),
        )

        point_model.save_as_nifti(
            field_mags_full,
            os.path.join(self.output_path, "VTA_solution.nii"),
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

    def h1_space(self, boundaries: List[str], is_complex: bool) -> ngsolve.H1:
        """Return a h1 space on the mesh.

        Parameters
        ----------
        boundaries : list of str
            List of boundary names.
        is_complex: bool
            Whether to use complex arithmetic

        Returns
        -------
        ngsolve.H1
        """
        dirichlet = "|".join(boundary for boundary in boundaries)
        return ngsolve.H1(
            mesh=self.mesh.ngsolvemesh,
            order=self._order,
            dirichlet=dirichlet,
            complex=is_complex,
            wb_withedges=False,
        )

    def number_space(self) -> ngsolve.comp.NumberSpace:
        """Return a number space on the mesh.

        Returns
        -------
        ngsolve.NumberSpace
            Space with only one single (global) DOF.

        TODO check if needed
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

    def get_scale_factor(self, freq_idx: int) -> float:
        """Scale solution by signal amplitude at a frequency given by index.

        Notes
        -----
        In voltage-controlled mode,
        only the amplitude of the Fourier coefficient is used.
        In current-controlled mode, TODO
        """
        scale_factor = 1.0
        if self.current_controlled:
            _logger.debug("Scale solution for current_controlled mode")
            if self.current_controlled and len(self.contacts.active) == 2:
                impedance = self.impedances[freq_idx]
                # use Ohm's law U = Z * I
                # and that the Fourier coefficient for the current is known
                amplitude = self.contacts.active[0].current
                # use positive current by construction
                if self.is_complex:
                    sign = np.sign(amplitude.real)
                else:
                    sign = np.sign(amplitude)
                amplitude *= sign
                scale_factor *= impedance * amplitude
        return scale_factor

    def prepare_current_controlled_mode(self) -> None:
        """Check contacts and assign voltages if needed."""
        if len(self.contacts.active) == 2:
            _logger.info("Overwrite voltage for current-controlled mode")
            for contact_idx, contact in enumerate(self.contacts.active):
                self.contacts[contact.name].voltage = float(contact_idx)
        else:
            if len(self.contacts.active) != 1:
                raise ValueError(
                    """In multicontact current-controlled mode,
                       currently only one active contact with fixed voltage can be used.
                       Its voltage has to be 0V (ground)."""
                )
            for contact in self.contacts.active:
                if not np.isclose(contact.voltage, 0):
                    raise ValueError(
                        """In multicontact current-controlled mode,
                           only ground voltage (0V) can be set on active contacts!"""
                    )

    def setup_timings_dict(self, export_vtk: bool, point_model: bool) -> dict:
        """Setup dictionary to save execution times estimate."""
        timings = {}
        timings["ComputeSolution"] = []
        if export_vtk:
            timings["VTKExport"] = []
        if point_model is not None:
            timings["CopyValues"] = []
        return timings

    def _store_solution_at_contacts(
        self, band_indices: Union[List, np.ndarray]
    ) -> None:
        """Save voltages / currents at given frequency band for all contacts."""
        if self.current_controlled:
            for contact_idx, contact in enumerate(self.contacts.active):
                for freq_idx in band_indices:
                    scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
                    self._free_stimulation_variable[freq_idx, contact_idx] = (
                        scale_factor * contact.voltage
                    )
                    self._stimulation_variable[freq_idx, contact_idx] = (
                        scale_factor * contact.current
                    )
        else:
            estimated_currents = self.estimate_currents()
            for contact_idx, contact in enumerate(self.contacts.active):
                for freq_idx in band_indices:
                    scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
                    self._free_stimulation_variable[freq_idx, contact_idx] = (
                        scale_factor * estimated_currents[contact.name]
                    )
                    self._stimulation_variable[freq_idx, contact_idx] = (
                        scale_factor * contact.voltage
                    )

    def _copy_values_for_time_domain(
        self,
        band_indices: Union[List, np.ndarray],
        tmp_array: np.ndarray,
        values: np.ndarray,
    ) -> None:
        """Copy values to time-domain vector."""
        for freq_idx in band_indices:
            scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
            tmp_array[freq_idx, :] = scale_factor * values

    def _has_sigma_changed(self, freq_idx, threshold=0.01) -> bool:
        """Check if conductivity has changed."""
        if self._sigma is None:
            return True
        else:
            max_error = 0.0
            for _material, model in self._conductivity_cf.dielectric_properties.items():
                if self.is_complex:
                    old_value = model.complex_conductivity(
                        (freq_idx - 1) * self.signal.base_frequency
                    )
                    new_value = model.complex_conductivity(
                        freq_idx * self.signal.base_frequency
                    )
                    error_real = np.abs(old_value.real - new_value.real)
                    # to catch zero-case
                    if not np.isclose(old_value.real, 0.0):
                        error_real /= old_value.real
                    error_imag = np.abs(old_value.imag - new_value.imag)
                    # to catch zero-case
                    if not np.isclose(old_value.imag, 0.0):
                        error_imag /= old_value.imag

                    error = np.maximum(error_real, error_imag)
                else:
                    old_value = model.conductivity(
                        (freq_idx - 1) * self.signal.base_frequency
                    )
                    new_value = model.conductivity(
                        freq_idx * self.signal.base_frequency
                    )
                    error = np.abs((old_value - new_value) / old_value)
                if error > max_error:
                    max_error = error
            if max_error > threshold:
                return True
            return False

    def _add_surface_impedance(self) -> bool:
        """Decide if surface impedance should be added to active contacts."""
        add_surface_impedance = False
        for contact in self.contacts.active:
            add_surface_impedance = not np.isclose(contact.surface_impedance, 0.0)
        return add_surface_impedance
