# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk
# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import logging
import os
import time
from abc import ABC, abstractmethod
from typing import Optional, Union

import ngsolve
import numpy as np
import pandas as pd

from ossdbs.fem.mesh import Mesh
from ossdbs.fem.solver import Solver
from ossdbs.model_geometry import Contacts, ModelGeometry
from ossdbs.point_analysis import Lattice, PointModel
from ossdbs.stimulation_signals import (
    FrequencyDomainSignal,
    get_indices_in_octave_band,
    get_octave_band_indices,
    get_timesteps,
    reconstruct_time_signals,
)
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
    output_path: str
        Path to store solution
    """

    def __init__(
        self,
        geometry: ModelGeometry,
        conductivity: ConductivityCF,
        solver: Solver,
        order: int,
        meshing_parameters: dict,
        output_path: str = "Results",
    ) -> None:
        self._solver = solver
        self._order = order

        self._model_geometry = geometry

        # contacts of electrode
        self._contacts = Contacts(geometry.contacts)

        _logger.debug(f"Assigned base contacts with properties:\n {self._contacts}")

        self._conductivity_cf = conductivity
        self._complex = conductivity.is_complex

        # to store impedances at all frequencies
        self._impedances = None
        # to store the voltage (current-controlled)
        # or current (voltage-controlled) stimulation
        self._free_stimulation_variable = None
        self._stimulation_variable = None
        self._floatings_potentials = None

        # set output path
        self.output_path = output_path

        # generate the mesh
        self._mesh = Mesh(self._model_geometry.geometry, self._order)
        if meshing_parameters["LoadMesh"]:
            self.mesh.load_mesh(meshing_parameters["LoadPath"])
        else:
            self.mesh.generate_mesh(meshing_parameters)

        if meshing_parameters["SaveMesh"]:
            self.mesh.save(meshing_parameters["SavePath"])

        # to save previous solution and do post-processing
        self._frequency = None
        self._sigma = None

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
        return

    @abstractmethod
    def update_space(self):
        """Update space (e.g., if mesh changes)."""

    # ruff: noqa: C901
    def run_full_analysis(
        self,
        frequency_domain_signal: FrequencyDomainSignal,
        compute_impedance: bool = False,
        export_vtk: bool = False,
        point_models: Optional[list[PointModel]] = None,
        activation_threshold: Optional[float] = None,
        out_of_core: bool = False,
        export_frequency: Optional[float] = None,
        adaptive_mesh_refinement_settings: Optional[dict] = None,
        material_mesh_refinement_steps: int = 0,
        truncation_time: Optional[float] = None,
    ) -> dict:
        """Run volume conductor model at all frequencies.

        Parameters
        ----------
        compute_impedance: bool
            If True, the impedance will be computed at each frequency.
        export_vtk: bool
            VTK export for visualization in ParaView
        point_models: list[PointModel]
            list of PointModel to extract solution for VTA / PAM
        activation_threshold: float
            If VTA is estimated by threshold, provide it here.
            Its unit must be V/m!
        out_of_core: bool
            Indicate whether point model shall be done out-of-core
        export_frequency: float
            Frequency at which the VTK file should be exported.
            Otherwise, median frequency is used.
        frequency_domain_signal: FrequencyDomainSignal
            Frequency-domain representation of stimulation signal
        adaptive_mesh_refinement_settings: dict
            Perform adaptive mesh refinement (only at first frequency)
        material_mesh_refinement_steps: int
            How often should elements with more than one material be refined
        truncation_time: float
            Time until which result will be written to hard drive

        Notes
        -----
        The volume conductor model is run at all frequencies
        and the time-domain signal is computed (if relevant).
        """
        self.signal = frequency_domain_signal

        # set voltages for two-contact current controlled stimulation
        # check if current-controlled mode is correctly prepared
        if self.current_controlled:
            self.prepare_current_controlled_mode()

        if point_models is None:
            # use empty list
            point_models = []
        timings = self.setup_timings_dict(export_vtk, point_models)

        dtype = float
        if self.is_complex:
            dtype = complex

        _do_AMR = False
        if adaptive_mesh_refinement_settings is not None:
            self._check_AMR_settings(adaptive_mesh_refinement_settings)
            if "Active" in adaptive_mesh_refinement_settings:
                if adaptive_mesh_refinement_settings["Active"]:
                    _do_AMR = True

        multisine_mode = np.all(np.isclose(self.signal.amplitudes, 1.0))

        # always compute impedance for CC with 2 contacts
        if self.current_controlled and len(self.contacts.active) == 2:
            _logger.info(
                "Set compute_impedance to True."
                "Impedance calculation is required for 2 contacts."
            )
            compute_impedance = True

        if self.signal.octave_band_approximation:
            frequency_indices = get_octave_band_indices(self.signal.frequencies)
            # add DC component
            if not np.isclose(self.signal.amplitudes[0], 0.0):
                frequency_indices = np.insert(frequency_indices, 0, 0)
        else:
            frequency_indices = np.arange(len(self.signal.frequencies))

        if export_frequency is None:
            middle_frequency_index = frequency_indices[int(len(frequency_indices) / 2)]
            self._export_frequency = self.signal.frequencies[middle_frequency_index]
            _logger.info(f"Set export frequency to {self._export_frequency}")
        else:
            self._export_frequency = export_frequency

        if not multisine_mode:
            self._free_stimulation_variable = np.zeros(
                shape=(len(self.signal.frequencies), len(self.contacts.active)),
                dtype=complex,
            )
            self._stimulation_variable = np.zeros(
                shape=(len(self.signal.frequencies), len(self.contacts.active)),
                dtype=complex,
            )
            self._floatings_potentials = np.zeros(
                shape=(len(self.signal.frequencies), len(self.contacts.floating)),
                dtype=complex,
            )

        if compute_impedance:
            self._impedances = np.ndarray(
                shape=(len(self.signal.frequencies)), dtype=dtype
            )

        _logger.info(
            f"Number of elements before material refinement:{self.mesh.ngsolvemesh.ne}"
        )
        self.refine_mesh_by_material(material_mesh_refinement_steps)
        _logger.info(
            f"Number of elements after material refinement:{self.mesh.ngsolvemesh.ne}"
        )

        for computing_idx, freq_idx in enumerate(frequency_indices):
            frequency = self.signal.frequencies[freq_idx]
            _logger.info(f"Computing at frequency: {frequency}")
            time_0 = time.time()

            # prepare storing at multiple frequencies
            if self.signal.octave_band_approximation:
                band_indices = get_indices_in_octave_band(
                    freq_idx, frequency_indices, len(self.signal.frequencies) - 1
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

                # refine only at first frequency
                if computing_idx == 0 and _do_AMR:
                    if not compute_impedance:
                        try:
                            impedance = self.compute_impedance()
                        except NotImplementedError:
                            _logger.info("Using power instead of impedance in AMR")
                            impedance = self.compute_power()
                    else:
                        impedance = self._impedances[freq_idx]
                    _logger.info(
                        "Number of elements before refinement:"
                        f"{self.mesh.ngsolvemesh.ne}"
                    )
                    error = 100.0
                    refinements = 0
                    tolerance = adaptive_mesh_refinement_settings["ErrorTolerance"]
                    max_iterations = adaptive_mesh_refinement_settings["MaxIterations"]
                    while error > tolerance and refinements < max_iterations:
                        self.adaptive_mesh_refinement()
                        # solve on refined mesh
                        self.compute_solution(frequency)
                        # check new impedance
                        try:
                            new_impedance = self.compute_impedance()
                        except NotImplementedError:
                            new_impedance = self.compute_power()
                        # error in percent
                        error = 100.0 * abs(impedance - new_impedance) / abs(impedance)
                        # update variables for loop
                        refinements += 1
                        impedance = new_impedance
                        _logger.info(
                            f"Adaptive refinement step {refinements}, "
                            f"error {error:.3f}%."
                        )
                    if compute_impedance:
                        # overwrite impedance values
                        self._impedances[band_indices] = impedance

                    _logger.info(
                        "Number of elements after refinement:"
                        f"{self.mesh.ngsolvemesh.ne}"
                    )
                    _logger.info(
                        "Adaptive mesh refinement converged after "
                        f"{refinements} refinement steps with an "
                        f"error in the impedance of {error:.3f}%"
                    )
            else:
                _logger.info(f"Skipped computation at {frequency} Hz")
                if compute_impedance:
                    # copy from previous frequency
                    impedance = self._impedances[computing_idx - 1]
                    self._impedances[band_indices] = impedance
            # scale factor: is one for VC and depends on impedance for other case
            self._scale_factor = self.get_scale_factor(freq_idx)
            _logger.debug(f"Scale factor: {self._scale_factor}")

            if not multisine_mode:
                self._store_solution_at_contacts(band_indices)
            if _logger.getEffectiveLevel() == logging.DEBUG:
                estimated_currents = self.estimate_currents()
                _logger.debug(
                    f"Estimated currents through contacts: {estimated_currents}"
                )

            time_1 = time.time()
            timings["ComputeSolution"].append(time_1 - time_0)
            time_0 = time_1

            _logger.info("Copy solution to point models")
            # initialise point models only after possible mesh change
            # AMR can change points in mesh
            if computing_idx == 0:
                for point_model in point_models:
                    point_model.output_path = self.output_path
                    point_model.prepare_VCM_specific_evaluation(
                        self.mesh, self.conductivity_cf
                    )
                    point_model.prepare_frequency_domain_data_structure(
                        self.signal.signal_length, out_of_core
                    )

            # copy solution to point models
            self._process_frequency_domain_solution(band_indices, point_models)
            time_1 = time.time()
            timings["CopyValues"].append(time_1 - time_0)
            time_0 = time_1

            # export frequency-domain solution at one frequency
            if np.isclose(frequency, self._export_frequency):
                _logger.info(f"Exporting at {self._export_frequency}")
                # save vtk
                if export_vtk:
                    self.vtk_export(freq_idx, multisine_mode)
                    time_1 = time.time()
                    timings["VTKExport"].append(time_1 - time_0)
                    time_0 = time_1
                # continue with frequency-domain exports
                self._frequency_domain_exports(
                    point_models, freq_idx, activation_threshold
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

        # export time domain solution if a proper signal has been passed
        _logger.info("Launching reconstruction of time domain")
        if len(self.signal.frequencies) > 1 and not multisine_mode:
            _logger.info("Reconstructing time-domain signal.")
            timesteps = get_timesteps(
                self.signal.cutoff_frequency,
                self.signal.base_frequency,
                self.signal.signal_length,
            )

            truncation_index = None
            if truncation_time is not None:
                timestep = timesteps[1] - timesteps[0]
                truncation_index = round(truncation_time / timestep)
            for point_model_idx, point_model in enumerate(point_models):
                # skip point models that are not considered in time domain
                if not point_model.time_domain_conversion:
                    continue

                (
                    potential_in_time,
                    Ex_in_time,
                    Ey_in_time,
                    Ez_in_time,
                ) = point_model.compute_solutions_in_time_domain(
                    self.signal.signal_length, convert_field=point_model.export_field
                )
                point_model.create_time_result(
                    timesteps,
                    potential_in_time,
                    Ex_in_time,
                    Ey_in_time,
                    Ez_in_time,
                    truncation_index=truncation_index,
                )

                time_1 = time.time()
                timings[f"ReconstructTimeSignals_PointModel_{point_model_idx}"] = (
                    time_1 - time_0
                )
                time_0 = time_1

        # close output-file
        # and write point model reports
        for point_model in point_models:
            point_model.close_output_file()
            try:
                point_model.export_point_model_information(
                    os.path.join(point_model.output_path, point_model.name + ".json")
                )
            except NotImplementedError:
                pass

        if len(self.signal.frequencies) > 1 and not multisine_mode:
            self.export_solution_at_contacts()
        self._save_report(timings)
        return timings

    def export_solution_at_contacts(self) -> None:
        """Time-domain solution export."""
        timesteps = get_timesteps(
            self.signal.cutoff_frequency,
            self.signal.base_frequency,
            self.signal.signal_length,
        )
        free_stimulation_variable_in_time = reconstruct_time_signals(
            self._free_stimulation_variable, self.signal.signal_length
        )
        stimulation_variable_in_time = reconstruct_time_signals(
            self._stimulation_variable, self.signal.signal_length
        )
        free_stimulation_variable_at_contact = {}
        free_stimulation_variable_at_contact["time"] = timesteps
        for contact_idx, contact in enumerate(self.contacts.active):
            free_stimulation_variable_at_contact[contact.name + "_free"] = (
                free_stimulation_variable_in_time[:, contact_idx]
            )
            free_stimulation_variable_at_contact[contact.name] = (
                stimulation_variable_in_time[:, contact_idx]
            )

        df = pd.DataFrame(free_stimulation_variable_at_contact)
        df.to_csv(
            os.path.join(self.output_path, "stimulation_in_time.csv"), index=False
        )

        floating_at_contact = {}
        floating_at_contact["time"] = timesteps
        floating_potentials_in_time = reconstruct_time_signals(
            self._floatings_potentials, self.signal.signal_length
        )
        for contact_idx, contact in enumerate(self.contacts.floating):
            floating_at_contact[contact.name] = floating_potentials_in_time[
                :, contact_idx
            ]
        df = pd.DataFrame(floating_at_contact)
        df.to_csv(os.path.join(self.output_path, "floating_in_time.csv"), index=False)

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
                active_contacts_grounded = np.isclose(voltages_active, 0.0)
                if len(np.where(active_contacts_grounded)[0]) != 1:
                    raise ValueError(
                        "In multipolar current-controlled mode, "
                        "only one active contact has to be grounded!"
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

    def compute_power(self) -> complex:
        """Compute power in domain."""
        mesh = self.mesh.ngsolvemesh
        # do not need to account for mm because of integration
        power = ngsolve.Integrate(
            ngsolve.Conj(self.electric_field) * self.current_density, mesh
        )
        return power

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
            power = self.compute_power()
            # TODO integrate surface impedance by thin layer
            voltage = 0
            for idx, contact in enumerate(self.contacts.active):
                voltage += (-1) ** idx * contact.voltage
            _logger.debug(f"Voltage drop for impedance: {voltage}")
            return voltage * np.conj(voltage) / power
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

    def vtk_export(self, freq_idx: int, multisine_mode: bool = False) -> None:
        """Export all relevant properties to VTK.

        Parameters
        ----------
        freq_idx: int
            Index of frequency
        multisine_mode: bool
            If rectangular pulse is used (multisine_mode = False)
        """
        self.export_solution_to_vtk(freq_idx, multisine_mode)
        self.export_conductivity_to_vtk()
        self.export_material_distribution_to_vtk()

    def export_solution_to_vtk(
        self, freq_idx: int, multisine_mode: bool = False
    ) -> None:
        """Export potential and field at frequency to VTK.

        Parameters
        ----------
        freq_idx: int
            Index of frequency
        multisine_mode: bool
            If rectangular pulse is used (multisine_mode = False)
        """
        ngmesh = self.mesh.ngsolvemesh
        # use standard solution with 1V voltage drop
        # unless we run multisine mode
        scale_factor = 1.0
        if multisine_mode:
            scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
        FieldSolution(
            scale_factor * self.potential, "potential", ngmesh, self.is_complex
        ).save(os.path.join(self.output_path, "potential"))

        FieldSolution(
            scale_factor * self.electric_field, "E_field", ngmesh, self.is_complex
        ).save(os.path.join(self.output_path, "E-field"))

    def export_conductivity_to_vtk(self) -> None:
        """Write conductivity to VTK file."""
        ngmesh = self.mesh.ngsolvemesh
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

    def export_material_distribution_to_vtk(self) -> None:
        """Write material distribution to VTK file."""
        ngmesh = self.mesh.ngsolvemesh
        FieldSolution(
            self.conductivity_cf.material_distribution(self.mesh),
            "material",
            ngmesh,
            False,
        ).save(os.path.join(self.output_path, "material"))

    def floating_values(self) -> dict:
        """Read out floating potentials."""
        floating_voltages = {}
        for contact in self.contacts.floating:
            floating_voltages[contact.name] = contact.voltage
        return floating_voltages

    def h1_space(self, boundaries: list[str], is_complex: bool) -> ngsolve.H1:
        """Return a h1 space on the mesh.

        Parameters
        ----------
        boundaries : list of str
            list of boundary names.
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
            mesh=self.mesh.ngsolvemesh,
            order=max(1, self._order - 1),
            complex=self.is_complex,
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
            ground_assigned = False
            for contact_idx, contact in enumerate(self.contacts.active):
                if np.isclose(self.contacts[contact.name].voltage, 0):
                    if ground_assigned:
                        raise ValueError(
                            "All active contacts have been grounded (voltage = 0V)."
                            "Choose only one ground."
                        )
                    ground_assigned = True
                else:
                    contact_voltage = float(contact_idx) + 1
                    self.contacts[contact.name].voltage = contact_voltage
        else:
            if len(self.contacts.active) != 1:
                raise ValueError(
                    "In multicontact current-controlled mode,"
                    "currently only one active contact with fixed voltage can be used."
                    "Its voltage has to be 0V (ground)."
                )
            for contact in self.contacts.active:
                if not np.isclose(contact.voltage, 0):
                    raise ValueError(
                        "In multicontact current-controlled mode,"
                        "only ground voltage (0V) can be set on active contacts!"
                    )

    def setup_timings_dict(
        self, export_vtk: bool, point_models: list[PointModel]
    ) -> dict:
        """Setup dictionary to save execution times estimate."""
        timings = {}
        timings["ComputeSolution"] = []
        # look at entire copying process only
        timings["CopyValues"] = []
        if export_vtk:
            timings["VTKExport"] = []
        for point_model_idx, _ in enumerate(point_models):
            timings[f"ReconstructTimeSignals_PointModel_{point_model_idx}"] = 0.0
        return timings

    def _store_solution_at_contacts(
        self, band_indices: Union[list, np.ndarray]
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
        for contact_idx, contact in enumerate(self.contacts.floating):
            for freq_idx in band_indices:
                scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
                self._floatings_potentials[freq_idx, contact_idx] = (
                    scale_factor * contact.voltage
                )

    def _copy_frequency_domain_solution(
        self,
        band_indices: Union[list, np.ndarray],
        point_model: PointModel,
        potentials: np.ndarray,
        fields: np.ndarray,
    ) -> None:
        """Copy values to time-domain vector."""
        for freq_idx in band_indices:
            scale_factor = self._scale_factor * self.signal.amplitudes[freq_idx]
            # cast scale_factor to complex
            # needed for export of h5py files in out-of-core mode
            if not isinstance(scale_factor, complex):
                scale_factor = complex(scale_factor)
            point_model.copy_frequency_domain_solution_from_vcm(
                freq_idx, scale_factor * potentials, scale_factor * fields
            )

    def threshold_frequency_domain_Efield(
        self, scale_factor: float, activation_threshold: float
    ) -> float:
        """Determine volume of E-field above threshold at current frequency."""
        field = scale_factor * self.electric_field
        # convert to V/m
        field_magnitude = 1000.0 * ngsolve.sqrt(ngsolve.InnerProduct(field, field))
        # subtract threshold from electric field,
        # all positive values are 1, negative values 0
        threshold_cf = ngsolve.IfPos(field_magnitude - activation_threshold, 1, 0)
        mesh = self.mesh.ngsolvemesh
        # Integrate to get volume
        return ngsolve.Integrate(threshold_cf, mesh=mesh)

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

    def _frequency_domain_exports(
        self,
        point_models: list,
        export_frequency_index: int,
        activation_threshold: float,
    ):
        """Export solution at desired frequency."""
        export_frequency = self.signal.frequencies[export_frequency_index]
        _logger.info(f"Exporting results at {export_frequency} Hz.")
        for point_model in point_models:
            _logger.info(f"Exporting for point model type {type(point_model)}.")
            point_model.export_potential_at_frequency(
                self._export_frequency, export_frequency_index
            )
            if point_model._export_field:
                _logger.info(
                    f"Exporting electric field at frequency {self._export_frequency}Hz."
                )
                point_model.export_field_at_frequency(
                    self._export_frequency,
                    export_frequency_index,
                    electrode=self.model_geometry.electrodes[0],
                    activation_threshold=activation_threshold,
                )
            if isinstance(point_model, Lattice):
                scale_factor = (
                    self._scale_factor * self.signal.amplitudes[export_frequency_index]
                )
                point_model.VTA_volume = self.threshold_frequency_domain_Efield(
                    scale_factor, activation_threshold
                )
                _logger.info(f"VTA volume is: {point_model.VTA_volume:.3f}")

    def _process_frequency_domain_solution(
        self, band_indices: Union[list, np.ndarray], point_models: PointModel
    ):
        """Copy results to points."""
        for point_model in point_models:
            potentials = self.evaluate_potential_at_points(point_model.lattice)
            fields = self.evaluate_field_at_points(point_model.lattice)

            # copy values for time-domain analysis
            self._copy_frequency_domain_solution(
                band_indices, point_model, potentials, fields
            )

    def adaptive_mesh_refinement(self):
        """Refine mesh adaptively."""
        flux = self.current_density
        hdiv_space = self.flux_space()
        flux_potential = ngsolve.GridFunction(space=hdiv_space)
        flux_potential.Set(coefficient=flux)
        difference = flux - flux_potential
        error = difference * ngsolve.Conj(difference)
        self.mesh.refine_by_error_cf(error)
        self.update_space()

    def _check_AMR_settings(self, adaptive_mesh_refinement_settings: dict) -> None:
        if not {"ErrorTolerance", "MaxIterations"}.issubset(
            adaptive_mesh_refinement_settings.keys()
        ):
            raise ValueError(
                "Need to specify ErrorTolerance and "
                "MaxIterations for adaptive mesh refinement"
            )

    def _save_report(self, timings: dict):
        """Save simulation run report to disk."""
        report = {}
        report["DOF"] = self._space.ndof
        report["Elements"] = self.mesh.n_elements
        report["Timings"] = timings

        with open(os.path.join(self.output_path, "VCM_report.json"), "w") as fp:
            json.dump(report, fp)

    def refine_mesh_by_material(self, material_mesh_refinement_steps: int) -> None:
        """Check materials and refine mesh if more than one material per element."""
        for _ in range(material_mesh_refinement_steps):
            self.mesh.refine_by_material_cf(
                self.conductivity_cf.material_distribution(self.mesh)
            )
