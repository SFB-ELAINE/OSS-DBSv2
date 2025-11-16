# Copyright 2023, 2024 Konstantin Butenko, Jan Philipp Payonk, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from dataclasses import asdict

import h5py
import numpy as np
import scipy

from ossdbs.electrodes.defaults import default_electrode_parameters


class LeadSettings:
    """Lead-DBS settings in OSS-DBS format.

    Parameters
    ----------
    mat_file_path : str
        lead-dbs settings created in ea_genvat_butenko.m

    Notes
    -----
    The lead-dbs settings are usually stored in oss-dbs_parameters.mat

    """

    def __init__(self, mat_file_path):
        # load .mat of different versions
        try:
            self._file = h5py.File(str(mat_file_path), "r")
            self._h5 = True
        except ValueError:
            print(
                "\n Please, save oss-dbs_parameters using"
                "'save(oss-dbs_parameters_path, 'settings', '-v7.3')'"
            )
            raise
            # TODO  Fix non-binary .mat import
            self._file = scipy.io.loadmat(mat_file_path)
            self._h5 = False

        self._settings = self._file["settings"]

        # Check that both electrodes are either CC or VC
        cur_ctrl_arr = self.get_cur_ctrl()
        if all(~np.isnan(cur_ctrl_arr)):
            if cur_ctrl_arr[0] != cur_ctrl_arr[1]:
                raise RuntimeError("Simultaneous use of VC and CC is not allowed!")

        # for now restrict to one electrode per simulation
        self.NUM_ELECS = 1

    def make_oss_settings(self, hemis_idx: int = 0, output_path: str = "") -> dict:
        """Convert lead settings into OSS-DBS parameters.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        output_path: str
            path to write file

        Returns
        -------
        dict
        """
        # Constants
        # OSS strings used for tissue types
        UNK, CSF, GM, WM, BLOOD = [
            "Unknown",
            "CSF",
            "Gray matter",
            "White matter",
            "Blood",
        ]

        # folders for the output
        HEMIS_OUTPUT_PATH = "Results"
        SIDES = ["rh", "lh"]
        side = SIDES[hemis_idx]

        (
            elec_dict,
            unit_directions,
            specs_array_length,
        ) = self.import_implantation_settings(hemis_idx)

        # define if current-controlled
        current_controlled = bool(self.get_cur_ctrl()[hemis_idx])

        # get electrode dictionary with stimulation parameters
        elec_dict, case_grounding, floating = self.import_stimulation_settings(
            hemis_idx, current_controlled, elec_dict
        )

        # add another dict if more than one electrode at a time
        elec_dicts = [elec_dict]

        # enforce case grounding for current-controlled
        if current_controlled:
            case_grounding = True
            # SHOULD BE REMOVED LATER, only V = 0 is enough
            pulses_sign_amplitude = self.get_phi_vec() * 0.001  # switch to A
            pulse_sign_amplitude = pulses_sign_amplitude[hemis_idx, :]
            grounded_current = -1 * np.round(
                np.sum(pulse_sign_amplitude[~np.isnan(pulse_sign_amplitude)]), 9
            )  # could be 0
        else:
            # otherwise not relevant, but set to 0.0 if non-active contacts present
            grounded_current = 0.0

        grid_center, grid_resolution = self.get_grid_parameters(
            elec_dict["Name"], hemis_idx, unit_directions, specs_array_length
        )

        # check the distance between first and last contact for brain approx. dimensions
        contact_locations = self.get_cntct_loc(hemis_idx)
        first_contact = np.array(
            [
                contact_locations[0][0],
                contact_locations[1][0],
                contact_locations[2][0],
            ]
        )
        last_contact = np.array(
            [
                contact_locations[0][-1],
                contact_locations[1][-1],
                contact_locations[2][-1],
            ]
        )
        actual_span = np.linalg.norm(last_contact - first_contact)

        # use a conductivity model from the GUI
        cond_model = self.get_cond_model()

        # MAKE THE DICTIONARY
        partial_dict = {
            "ModelSide": 0,  # hardcoded for now, always keep to 0
            "BrainRegion": {
                # center at the head marker
                "Center": {
                    "x[mm]": self.get_imp_coord()[hemis_idx, 0],
                    "y[mm]": self.get_imp_coord()[hemis_idx, 1],
                    "z[mm]": self.get_imp_coord()[hemis_idx, 2],
                },
                # define brain approximation according to the electrode directionality
                "Dimension": {
                    "x[mm]": 50.0
                    + np.abs(unit_directions[hemis_idx, 0]) * actual_span * 2.0,
                    "y[mm]": 50.0
                    + np.abs(unit_directions[hemis_idx, 1]) * actual_span * 2.0,
                    "z[mm]": 50.0
                    + np.abs(unit_directions[hemis_idx, 2]) * actual_span * 2.0,
                },
            },
            "Electrodes": elec_dicts,
            "Surfaces": [
                {
                    "Name": "BrainSurface",
                    "Active": bool(case_grounding),
                    "Current[A]": grounded_current,
                    "Voltage[V]": 0.0,
                }
            ],
            "MaterialDistribution": {
                "MRIPath": self.get_mri_name(),
                "MRIMapping": {
                    # Make UNK map to one lower than the indices provided in
                    # the lead settings
                    UNK: int(
                        min([self.get_csf_idx(), self.get_gm_idx(), self.get_wm_idx()])
                        - 1
                    ),
                    CSF: self.get_csf_idx(),
                    WM: self.get_wm_idx(),
                    GM: self.get_gm_idx(),
                    # Make BLOOD map to one higher than the indices provided
                    # in the lead settings
                    BLOOD: int(
                        max([self.get_csf_idx(), self.get_gm_idx(), self.get_wm_idx()])
                        + 1
                    ),
                },
                "DiffusionTensorActive": len(self.get_dti_name()) > 0,
                "DTIPath": self.get_dti_name(),
            },
            "DielectricModel": {"Type": cond_model, "CustomParameters": None},
            "PointModel": {
                "Lattice": {
                    "Active": not (
                        bool(self.get_calc_axon_act())
                    ),  # for now, disable lattice when using PAM
                    "Center": {
                        # center at the middle of the electrode array
                        "x[mm]": grid_center[0],
                        "y[mm]": grid_center[1],
                        "z[mm]": grid_center[2],
                    },
                    "Shape": {"x": 71, "y": 71, "z": 71},
                    "Direction": {
                        "x[mm]": 0,
                        "y[mm]": 0,
                        "z[mm]": 1,
                    },
                    "PointDistance[mm]": grid_resolution,
                    "CollapseVTA": bool(self.remove_electrode()),
                },
                "Pathway": {
                    "Active": bool(self.get_calc_axon_act()),
                    "FileName": os.path.join(output_path, "Allocated_axons.h5"),
                },
            },
            "StimulationSignal": {"CurrentControlled": current_controlled},
            "CalcAxonActivation": bool(self.get_calc_axon_act()),
            "ActivationThresholdVTA[V-per-m]": float(
                self.get_act_thresh_vta()[hemis_idx]
            ),
            "OutputPath": os.path.join(output_path, HEMIS_OUTPUT_PATH),
            "FailFlag": side,
            "TemplateSpace": self.get_est_in_temp(),
            "Solver": {},
            "OutOfCore": bool(self.get_out_of_core()),
            "StimSets": {
                "Active": bool(self.get_stim_set_mode()),
                "StimSetsFile": os.path.join(
                    output_path, f"Current_protocols_{hemis_idx}.csv"
                ),
            },
        }

        # the constant model from Lead-DBS also assumes homogeneity
        # by default, gray matter parameters at 10 kHz, but can be customized
        if cond_model == "Constant":
            partial_dict["DielectricModel"]["CustomParameters"] = {
                "Gray matter": {
                    "permittivity": 2.22e4,
                    "conductivity": 0.2,
                },
                "White matter": {
                    "permittivity": 2.22e4,
                    "conductivity": 0.2,
                },
                "CSF": {
                    "permittivity": 2.22e4,
                    "conductivity": 0.2,
                },
                "Blood": {
                    "permittivity": 2.22e4,
                    "conductivity": 0.2,
                },
                "Unknown": {
                    "permittivity": 2.22e4,
                    "conductivity": 0.2,
                },
            }
        else:
            cond_model = self.get_cond_model()

        # use actual signal parameters for PAM
        if self.get_calc_axon_act():
            partial_dict = self.add_stimsignal_params(partial_dict, hemis_idx)
            # add path to the pathway parameter file
            partial_dict["PathwayFile"] = os.path.join(
                output_path, self.get_pathway_params_file()
            )

        # do not use h1amg as coarsetype preconditioner
        # if floating potentials are involved
        # set also to floating if multicontact current-controlled
        if not floating:
            if current_controlled:
                floating = True
        if floating:
            partial_dict["Solver"]["Preconditioner"] = "local"
            # increase number of iterations for FFEM
            if partial_dict["CalcAxonActivation"]:
                partial_dict["Solver"]["MaximumSteps"] = 2000

        return partial_dict

    def get_num_elecs(self):
        """Number of electrodes."""
        return self.NUM_ELECS

    def get_pat_fold(self):
        """Patient folder."""
        return self._get_str("Patient_folder")

    def get_est_in_temp(self):
        """TODO description."""
        return bool(self._get_num("Estimate_In_Template"))

    def get_mri_name(self):
        """Path to MRI file."""
        mri_file = self._get_str("MRI_data_name")
        return os.path.abspath(mri_file)

    def get_dti_name(self):
        """Path to DTI file."""
        return self._get_str("DTI_data_name")

    def get_cond_model(self):
        """Conductivity Model."""
        return self._get_str("cond_model")

    def get_gm_idx(self):
        """Gray matter index in MRI."""
        return int(self._get_num("GM_index"))

    def get_wm_idx(self):
        """White matter index in MRI."""
        return int(self._get_num("WM_index"))

    def get_csf_idx(self):
        """CSF index in MRI."""
        return int(self._get_num("CSF_index"))

    def get_def_mat(self):
        """Default material."""
        return self._get_str("default_material")

    def get_elec_type(self):
        """Electrode type."""
        return self._get_str("Electrode_type")

    # TODO check if calculation matches Lead-DBS (Konstantin)
    def get_specs_array_length(self, oss_elec_name):
        """TODO description."""
        elec_params = default_electrode_parameters[oss_elec_name]
        return elec_params.get_distance_l1_l4()

    def get_stretch_factor(self, specs_array_length, hemis_idx=0) -> float:
        """Compute stretch/squeeze factor between the first and the last
        contact. Relevant for normative space computations since the
        electrodes are non-linearly warped.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        specs_array_length: float
            distance between the first and the last contact

        Returns
        -------
        stretch_factor: float
            describes stretching for the electrode
        """
        C1_coords = self.get_imp_coord()[hemis_idx, :]
        C_last_coords = self.get_sec_coord()[hemis_idx, :]
        el_array_length = np.sqrt(
            (C1_coords[0] - C_last_coords[0]) ** 2
            + (C1_coords[1] - C_last_coords[1]) ** 2
            + (C1_coords[2] - C_last_coords[2]) ** 2
        )
        return el_array_length / specs_array_length

    def get_cntct_loc(self, hemis_idx):
        """Contact location."""
        contactCoords = np.asarray(
            self._file[self._settings["contactLocation"][hemis_idx, 0]][:, :]
        )
        return contactCoords

    # Used to re-compute rot_z
    def get_y_mark_nat(self):
        """Native y-marker."""
        return self._get_arr("yMarkerNative")

    def get_y_mark_mni(self):
        """MNI y-marker."""
        return self._get_arr("yMarkerMNI")

    def get_head_nat(self):
        """Native head."""
        return self._get_arr("headNative")

    def get_head_mni(self):
        """MNI head."""
        return self._get_arr("headMNI")

    def get_imp_coord(self):
        """Implantation coordinate."""
        return self._get_arr("Implantation_coordinate")

    def get_sec_coord(self):
        """Second coordinate."""
        return self._get_arr("Second_coordinate")

    def get_stim_center(self):
        """Stimulation center."""
        return self._get_arr("stim_center")

    def get_rot_z(self, index_side: int):
        """TODO description.

        Parameters
        ----------
        index_side: int
            Side of brain

        Notes
        -----
         Always recalculated from the other settings
        *IMPORTANT*: it is actually not native but scrf!
        """
        if self.get_est_in_temp():
            head_MNI = self.get_head_mni()[index_side, :]
            y = self.get_y_mark_mni()[index_side, :] - head_MNI
        else:
            head_nat = self.get_head_nat()[index_side, :]
            y = self.get_y_mark_nat()[index_side, :] - head_nat
        y_postop = y / np.linalg.norm(y)
        phi = np.arctan2(-y_postop[0], y_postop[1])
        return phi * 180.0 / np.pi

    def get_stim_set_mode(self):
        """Stimulation mode."""
        return self._get_num("stimSetMode")

    def get_cur_ctrl(self):
        """Current-controlled."""
        return self._get_arr("current_control").T[0]

    def remove_electrode(self):
        """Remove electrode (collapse VTA?)."""
        return self._get_num("removeElectrode")

    def get_neuron_model(self):
        """Neuron model."""
        return self._get_str("neuronModel")

    def get_signal_type(self):
        """Signal type."""
        return self._get_str("signalType")

    def get_pulse_width(self):
        """Pulse width."""
        return self._get_arr("pulseWidth")

    def check_biphasic(self):
        """Biphasic pulse."""
        return self._get_num("biphasic")

    def do_adaptive_ref(self):
        """Adaptive mesh refinement."""
        return self._get_num("AdaptiveRef")

    def get_encapsulation_type(self):
        """Encapsulation type."""
        return self._get_str("encapsulationType")

    def get_act_thresh_vta(self):
        """Activation threshold for VTA."""
        return self._get_arr("Activation_threshold_VTA")

    def get_phi_vec(self):
        """TODO description."""
        return self._get_arr("Phi_vector")

    def get_case_grnd(self):
        """Case grounding."""
        return self._get_arr("Case_grounding")

    def get_calc_axon_act(self):
        """Calculate axon activation."""
        return self._get_num("calcAxonActivation")

    def get_out_of_core(self):
        """Check if intermediate solution is unloaded."""
        return self._get_num("outOfCore")

    def get_connectome(self):
        """TODO description."""
        return self._get_str("connectome")

    def get_axon_len(self):
        """Axon length."""
        return self._get_arr("axonLength")

    def get_fib_diam(self):
        """Fibre diameter."""
        return self._get_arr("fiberDiameter")

    def get_conectome_path(self):
        """Connectome path."""
        return self._get_str("connectomePath")

    def get_pathway_params_file(self):
        """File with the pathway parameters file."""
        return self._get_str("pathwayParameterFile")

    def get_connectome_tract_names(self):
        """Get tract names in the connectome.

        Returns
        -------
        list
        """
        ref_arr = self._settings["connectomeTractNames"][0]

        # Initialize the tract names array
        tract_names = np.empty(len(ref_arr), dtype=object)

        # Convert the references to name strings and return the strings
        for i in range(len(ref_arr)):
            entry = self._file[ref_arr[i]][:]
            ascii_codes = np.ndarray(entry.shape[0])
            for j in range(entry.shape[0]):
                ascii_codes[j] = entry[j][0]
            tract_names[i] = "".join(np.vectorize(chr)(ascii_codes.astype(int)))
        return tract_names

    def get_inter_mode(self):
        """Interactive mode."""
        return self._get_num("interactiveMode")

    def add_stimsignal_params(self, partial_dict: dict, hemi_idx: int):
        """Add stimulation signal parameters."""
        partial_dict["StimulationSignal"]["Type"] = self.get_signal_type()
        if partial_dict["StimulationSignal"]["Type"] == "Train":
            partial_dict["StimulationSignal"]["Type"] = "Rectangle"
        partial_dict["StimulationSignal"]["PulseWidth[us]"] = float(
            self.get_pulse_width()[hemi_idx]
        )

        if self.check_biphasic():
            partial_dict["StimulationSignal"]["CounterPulseWidth[us]"] = float(
                self.get_pulse_width()[hemi_idx]
            )

        # hardwired for now
        partial_dict["StimulationSignal"]["Frequency[Hz]"] = 130.0
        partial_dict["StimulationSignal"]["SpectrumMode"] = "OctaveBand"
        partial_dict["StimulationSignal"]["CutoffFrequency"] = (
            250000.0  # 2 us time step
        )
        partial_dict["StimulationSignal"]["PulseTopWidth[us]"] = 0.0
        partial_dict["StimulationSignal"]["InterPulseWidth[us]"] = 0.0

        return partial_dict

    def stretch_electrode(self, oss_electrode_name: str, hemi_idx: int):
        """Stretch electrode geometry.

        Parameters
        ----------
        oss_electrode_name: str
            electrode name in OSS-DBS format
        hemi_idx: int
            Index of brain side
        """
        stretch_list = ["tip_length", "contact_length", "contact_spacing"]
        specs_array_length = self.get_specs_array_length(oss_electrode_name)
        stretch_factor = self.get_stretch_factor(specs_array_length, hemi_idx)
        if abs(stretch_factor - 1.0) < 0.01:  # 1% tolerance
            stretch_factor = 1.0
        default_parameters = default_electrode_parameters[oss_electrode_name]
        stretched_parameters = {}
        for parameter, value in zip(
            asdict(default_parameters), asdict(default_parameters).values()
        ):
            if parameter in stretch_list:
                stretched_parameters[parameter] = value * stretch_factor
            else:
                stretched_parameters[parameter] = value
        return stretched_parameters

    def get_tip_position(self, oss_elec_name: str, hemi_idx: int):
        """Get tip, implantation trajectory from head
        (Implantation_coordinate) and tail (Second_coordinate),
        and length of the contact span.

        Parameters
        ----------
        oss_elec_name: str
            electrode name in OSS-DBS format
        hemi_idx: int
            Index of brain side

        Returns
        -------
        numpy.ndarray, numpy.ndarray, float
        """
        elec_params = default_electrode_parameters[oss_elec_name]
        imp_coords = np.array(self.get_imp_coord())
        sec_coords = np.array(self.get_sec_coord())
        directions = sec_coords - imp_coords
        unit_directions = directions / np.linalg.norm(directions, axis=1)[:, np.newaxis]
        specs_array_length = self.get_specs_array_length(oss_elec_name)
        stretch_factor = self.get_stretch_factor(specs_array_length, hemi_idx)

        if abs(stretch_factor - 1.0) < 0.01:  # 1% tolerance
            stretch_factor = 1.0
        offset = elec_params.get_center_first_contact() * stretch_factor
        tip_position = imp_coords - offset * unit_directions

        return unit_directions, tip_position, specs_array_length

    def get_grid_parameters(
        self, electrode_name, hemis_idx, unit_directions, specs_array_length
    ):
        """Center lattice w.r.t. etimated stimulation volume and set resolution.

        Parameters
        ----------
        electrode_name: str
            OSS-DBS notation
        unit_directions: numpy.ndarray
            implantation trajectory
        specs_array_length: float
            length of the contact span
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)

        Returns
        -------
        grid_center: numpy.ndarray, center of the lattice model
        grid_resolution: float, resolution of the lattice model
        """
        # get grid center for lattice / voxel lattice model
        if np.any(np.isnan(self.get_stim_center()[hemis_idx, :])):
            self.grid_center = (
                self.get_imp_coord()[hemis_idx, :]
                + unit_directions[hemis_idx, :] * specs_array_length / 2
            )
        else:
            grid_center = self.get_stim_center()[hemis_idx, :]

        # set resolution
        # coarser resolution for large span electrodes
        # and large amplitudes (>5 mA or 5 V)
        phi_vector = self.get_phi_vec()[hemis_idx, :]
        if (
            electrode_name == "BostonScientificVercise"
            or electrode_name == "BostonScientificVerciseCustom"
            or electrode_name == "SceneRay1202"
            or electrode_name == "SceneRay1202Custom"
            or electrode_name == "SceneRay1212"
            or electrode_name == "SceneRay1212Custom"
            or electrode_name == "SceneRay1242"
            or electrode_name == "SceneRay1242Custom"
            or np.max(np.abs(phi_vector[~np.isnan(phi_vector)])) > 5.0
        ):
            grid_resolution = 0.4
        else:
            grid_resolution = 0.33

        return grid_center, grid_resolution

    def import_implantation_settings(self, hemis_idx, elec_dict=None):
        """Convert Lead-DBS implantation settings to OSS-DBS parameters.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        elec_dict: dict, Optional
            electrode dictionary to create/update

        Returns
        -------
        elec_dict: dict
        case_grounding: bool
        """
        # Convert the electrode name to OSS-DBS format
        electrode_name = self.get_elec_type()
        # Cartesia X not available is OSS-DBS, SNEX not available in Lead
        electrode_names = {
            "Abbott Directed 6172 (short)": "AbbottStJudeDirected6172",
            "Abbott Directed 6173 (long)": "AbbottStJudeDirected6173",
            "Abbott ActiveTip (6146-6149)": "AbbottStJudeActiveTip6146_6149",
            "Abbott ActiveTip (6142-6145)": "AbbottStJudeActiveTip6142_6145",
            "St. Jude ActiveTip (6146-6149)": "AbbottStJudeActiveTip6146_6149",
            "St. Jude ActiveTip (6142-6145)": "AbbottStJudeActiveTip6142_6145",
            "St. Jude Directed 6180": "AbbottStJudeDirected6172",
            "St. Jude Directed 6172 (short)": "AbbottStJudeDirected6172",
            "St. Jude Directed 6173 (long)": "AbbottStJudeDirected6173",
            "Boston Scientific Vercise": "BostonScientificVercise",
            "Boston Scientific Vercise Directed": "BostonScientificVerciseDirected",
            "Boston Scientific Vercise Cartesia HX": "BostonScientificCartesiaHX",
            "Boston Scientific Vercise Cartesia X": "BostonScientificCartesiaX",
            "ELAINE Rat Electrode": "MicroProbesRodentElectrode",
            "Medtronic 3387": "Medtronic3387",
            "Medtronic 3389": "Medtronic3389",
            "Medtronic 3391": "Medtronic3391",
            "SceneRay SR1210": "Medtronic3387",
            "SceneRay SR1200": "Medtronic3389",
            "SceneRay SR1242": "SceneRay1242",
            "SceneRay SR1202": "SceneRay1202",
            "SceneRay SR1212": "SceneRay1212",
            "SceneRay SR1211": "SceneRay1211",
            "Medtronic B33005": "MedtronicSenSightB33005",
            "Medtronic B33015": "MedtronicSenSightB33015",
            "PINS Medical L301": "PINSMedicalL301",
            "PINS Medical L302": "PINSMedicalL302",
            "PINS Medical L303": "PINSMedicalL303",
            "NeuroPace DL-344-3.5": "NeuroPaceDL344_3_5",
            "NeuroPace DL-344-10": "NeuroPaceDL344_10",
            "DIXI D08-05AM": "DixiSEEG5",
            "DIXI D08-08AM": "DixiSEEG8",
            "DIXI D08-10AM": "DixiSEEG10",
            "DIXI D08-12AM": "DixiSEEG12",
            "DIXI D08-15AM": "DixiSEEG15",
            "DIXI D08-18AM": "DixiSEEG18",
            "PMT 2102-08-091": "PMTsEEG2102_08",
            "PMT 2102-10-091": "PMTsEEG2102_10",
            "PMT 2102-12-091": "PMTsEEG2102_12",
            "PMT 2102-14-091": "PMTsEEG2102_14",
            "PMT 2102-16-091": "PMTsEEG2102_16",
            "PMT 2102-08-094": "PMTsEEG2102_08",  # same, just with permanent stylet
            "PMT 2102-10-094": "PMTsEEG2102_10",
            "PMT 2102-12-094": "PMTsEEG2102_12",
            "PMT 2102-14-094": "PMTsEEG2102_14",
            "PMT 2102-16-094": "PMTsEEG2102_16",
        }

        for lead in electrode_names.keys():
            if lead == electrode_name:
                electrode_name = electrode_names[lead]

        # Check that oss electrode name is valid
        if electrode_name not in default_electrode_parameters.keys():
            raise Exception(electrode_name + " is not a recognized electrode type")

        stretched_parameters = self.stretch_electrode(electrode_name, hemis_idx)

        # get tip position from the head and tail markers
        unit_directions, tip_pos, specs_array_length = self.get_tip_position(
            electrode_name, hemis_idx
        )

        elec_dict_imp = {
            # Assuming both electrodes are the same type
            "Name": electrode_name + "Custom",
            "CustomParameters": stretched_parameters,
            "Rotation[Degrees]": self.get_rot_z(hemis_idx),
            "Direction": {
                "x[mm]": unit_directions[hemis_idx, 0],
                "y[mm]": unit_directions[hemis_idx, 1],
                "z[mm]": unit_directions[hemis_idx, 2],
            },
            "TipPosition": {
                "x[mm]": tip_pos[hemis_idx, 0],
                "y[mm]": tip_pos[hemis_idx, 1],
                "z[mm]": tip_pos[hemis_idx, 2],
            },
            "EncapsulationLayer": {
                "Thickness[mm]": 0.1,
                "Material": self.get_encapsulation_type(),
                "DielectricModel": self.get_cond_model(),
                "DielectricParameters": None,
                "MaxMeshSize": 0.1,
            },
        }
        if elec_dict is None:
            elec_dict = elec_dict_imp
        else:
            elec_dict.update(elec_dict_imp)

        if elec_dict_imp["EncapsulationLayer"]["Material"] == "None":
            elec_dict_imp["EncapsulationLayer"]["Material"] = "Gray matter"
            elec_dict_imp["EncapsulationLayer"]["Thickness[mm]"] = 0.0

        return elec_dict, unit_directions, specs_array_length

    def import_stimulation_settings(self, hemis_idx, current_controlled, elec_dict):
        """Convert Lead-DBS stim settings to OSS-DBS parameters,
        update electrode dictionary.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        current_controlled: bool
            Current-controlled stimulation
        elec_dict: dict
            electrode dictionary to update

        Returns
        -------
        elec_dict: dict
        case_grounding: bool
        """
        # store if floating or not
        floating = False
        # stimulation vector over electrode contacts (negative - cathode)
        pulse_amps = self.get_phi_vec()

        # check if the stimulation was defined for this side
        if np.isnan(self.get_cur_ctrl()[hemis_idx]):
            raise RuntimeError("No stimulation defined for this side")

        # check stimulation mode
        if self.get_cur_ctrl()[hemis_idx]:
            # Lead-DBS uses mA as the input
            pulse_amps = pulse_amps * 0.001
        else:
            # Fix of VC random grounding bug for Lead-DBS stim settings
            pulse_amps[pulse_amps == 0] = float("nan")

        # make list of dictionaries for the electrode settings
        # for now use one electrode at a time

        for index_side in [hemis_idx]:
            pulse_amp = pulse_amps[index_side, :]

            if self.get_cur_ctrl()[index_side]:
                # for CC, check if currents sum up to 0.0.
                # If not, enable case grounding
                case_grounding = np.sum(pulse_amp[~np.isnan(pulse_amp)]) != 0
            else:
                # for VC, case grounding is defined explicitly
                case_grounding = bool(self.get_case_grnd()[index_side])

                # shift all voltages if bipolar case
                # to have 0V and cathodes (as in the stimulators)
                # but don't shift if purely anodic stim
                cathodic_case = case_grounding and np.nanmin(pulse_amp) >= 0.0
                if not cathodic_case and np.nanmax(pulse_amp) > 0.0:
                    pulse_amp[:] = pulse_amp[:] - np.nanmax(pulse_amp)

            # cntct_dicts is a list of the contacts that will go in the json
            # for this electrode

            cntct_dicts = np.empty(len(pulse_amp), dtype=object)
            cntcts_made = 0

            # get edge size from lead diameter
            # used in mesh refinement

            if "Name" not in elec_dict:
                raise KeyError("Need to provide name of electrode")
            electrode_name = elec_dict["Name"].replace("Custom", "")
            lead_diameter = default_electrode_parameters[electrode_name].lead_diameter
            perimeter = np.pi * lead_diameter
            edge_size = perimeter / 50.0

            for i in range(len(pulse_amp)):
                # all (truly) non-active contacts are floating with 0A
                if np.isnan(pulse_amp[i]):
                    floating = True
                    cntct_dicts[cntcts_made] = {
                        # Assuming one-indexed contact ids
                        "Contact_ID": i + 1,
                        "Active": False,
                        "Current[A]": 0.0,
                        "Voltage[V]": 0.0,
                        "Floating": True,
                    }
                else:
                    # for current-controlled, we have a pseudo non-active contact
                    if current_controlled:
                        cntct_dicts[cntcts_made] = {
                            # Assuming one-indexed contact ids
                            "Contact_ID": i + 1,
                            "Active": False,
                            "Current[A]": pulse_amp[i],
                            "Voltage[V]": 0.0,
                            "Floating": True,
                            "MaxMeshSizeEdge": edge_size,
                        }
                    else:
                        cntct_dicts[cntcts_made] = {
                            # Assuming one-indexed contact ids
                            "Contact_ID": i + 1,
                            "Active": True,
                            "Current[A]": 0.0,
                            "Voltage[V]": pulse_amp[i],
                            "Floating": False,
                            "MaxMeshSizeEdge": edge_size,
                        }

                cntcts_made += 1

        elec_dict["Contacts"] = cntct_dicts.tolist()

        return elec_dict, case_grounding, floating

    # Private fxns

    def _is_h5(self):
        return self._h5

    # Worth considering making different subclasses for different file types?
    def _get_num(self, field_name):
        if self._is_h5():
            return self._settings[field_name][0][0]
        else:
            return self._settings[field_name][0][0][0][0]

    def _get_str(self, field_name):
        entry = ""
        if self._is_h5():
            ascii_codes = self._settings[field_name][:, 0]
            entry = "".join(np.vectorize(chr)(ascii_codes))
        else:
            entry = self._settings[field_name][0][0][0]
        if entry == "no dti":
            entry = ""
        return entry

    def _get_arr(self, field_name):
        if self._is_h5():
            return self._settings[field_name][:, :].T
        else:
            if len(self._settings[field_name][0][0]) == 1:
                return self._settings[field_name][0][0][0]
            else:
                return self._settings[field_name][0][0].astype(float)
