# Copyright 2024 Konstantin Butenko, Julius Zimmermann
# Copyright 2017 Christian Schmidt
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import logging
import math
import os
from abc import ABC, abstractmethod
from typing import Optional

import h5py
import numpy as np
import scipy

from .axon_default_MRG2002 import get_axon_parameters_template
from .utilities import (
    convert_fibers_to_streamlines,
    normalized,
    place_axons_on_streamlines,
    resample_fibers_to_Ranviers,
)

_logger = logging.getLogger(__name__)


class AxonMorphology(ABC):
    """Axon morphology class."""

    def __init__(self, downsampled=False):
        self.downsampled = downsampled

        # defaults
        self._n_Ranvier = None

    @property
    def n_segments(self):
        """Number of segments."""
        return (self.n_Ranvier - 1) * self.n_comp + 1

    @property
    def fiber_diam(self):
        """Fiber diameter."""
        return self._fiber_diam

    @fiber_diam.setter
    def fiber_diam(self, value):
        if value < 0:
            raise ValueError("Fiber diameter has to be greater than zero")
        self._fiber_diam = value

    @property
    def axon_length(self):
        """Length of axon."""
        return self._axon_length

    @axon_length.setter
    def axon_length(self, value):
        self._axon_length = value
        expected_n_Ranvier = int(value / self.node_step)
        if not self.n_Ranvier == expected_n_Ranvier:
            # Updating number of Ranviers to match axon length
            self.n_Ranvier = expected_n_Ranvier

    @property
    def n_Ranvier(self):
        """Number of nodes of Ranvier compartments."""
        return self._n_Ranvier

    @n_Ranvier.setter
    def n_Ranvier(self, value: int):
        # must be an odd number!
        if value % 2 == 0:
            value -= 1
            # update axon length
            self.axon_length = value * self.node_step
        self._n_Ranvier = value

    @property
    def downsampled(self):
        """Full or downsampled model."""
        return self._downsampled

    @downsampled.setter
    def downsampled(self, value):
        self._downsampled = value

    @abstractmethod
    def update_axon_morphology(self, axon_diam, axon_length=None, n_Ranvier=None):
        """Axon morphology for specific model."""

    @abstractmethod
    def get_local_compartment_coords(self, axon_morphology):
        """Get 1-D coordinates of internodal compartments relative to the node at 0.0.

        Returns
        -------
        numpy.ndarray

        """

    @property
    def n_comp(self):
        """Number of compartments."""
        return self._n_comp

    @property
    @abstractmethod
    def node_step(self):
        """Node step."""


class AxonMorphologyMRG2002(AxonMorphology):
    """Axon morphology for MRG2002 model."""

    @AxonMorphology.n_comp.setter
    def n_comp(self, value):
        """Number of compartments."""
        self._n_comp = value

    @property
    def node_step(self):
        """Node step."""
        return self._node_step

    @node_step.setter
    def node_step(self, value):
        self._node_step = value

    @property
    def ranvier_length(self):
        """Length of Ranvier compartment."""
        return self._ranvier_length

    @ranvier_length.setter
    def ranvier_length(self, value):
        self._ranvier_length = value

    @property
    def para1_length(self):
        """Length of the 1st paranodal compartments."""
        return self._para1_length

    @para1_length.setter
    def para1_length(self, value):
        self._para1_length = value

    @property
    def para2_length(self):
        """Length of the 2nd paranodal compartments."""
        return self._para2_length

    @para2_length.setter
    def para2_length(self, value):
        self._para2_length = value

    @property
    def ranvier_nodes(self):
        """Number of node of Ranvier compartments."""
        return self._ranvier_nodes

    @ranvier_nodes.setter
    def ranvier_nodes(self, value):
        self._ranvier_nodes = int(value)

    @property
    def inter_nodes(self):
        """Number of internodal compartments."""
        return self._inter_nodes

    @inter_nodes.setter
    def inter_nodes(self, value):
        self._inter_nodes = int(value)

    @property
    def para1_nodes(self):
        """Number of first paranodal compartments."""
        return self._para1_nodes

    @para1_nodes.setter
    def para1_nodes(self, value):
        self._para1_nodes = int(value)

    @property
    def para2_nodes(self):
        """Number of 2nd paranodal compartments."""
        return self._para2_nodes

    @para2_nodes.setter
    def para2_nodes(self, value):
        self._para2_nodes = int(value)

    # TODO same as para1_nodes ?
    @property
    def n_para1(self):
        """Number of all 1st paranodal compartments."""
        return self._n_para1

    @n_para1.setter
    def n_para1(self, value):
        self._n_para1 = value

    @property
    def n_para2(self):
        """Number of all 2nd paranodal compartments."""
        return self._n_para2

    @n_para2.setter
    def n_para2(self, value):
        self._n_para2 = value

    @AxonMorphology.n_Ranvier.setter
    def n_Ranvier(self, value):
        """Nodes of Ranvier."""
        # must be an odd number!
        if value % 2 == 0:
            value -= 1
            # update axon length
            self.axon_length = value * self.node_step
        self._n_Ranvier = value
        # MRG2002 specific code
        self.n_para1 = (self.para1_nodes * self.n_Ranvier - 1) / (21 - 1)
        self.n_para2 = (self.para2_nodes * self.n_Ranvier - 1) / (21 - 1)

    def get_n_comp(self, downsampled):
        """Get number of compartments depending on sampling."""
        if downsampled:
            if self.fiber_diam >= 5.7:
                # node -- -- internodal -- -- -- -- internodal -- -- node
                n_comp = 3
            else:
                # node -- -- -- internodal -- -- -- node
                n_comp = 2
        else:
            n_comp = int(
                (
                    self.ranvier_nodes
                    - 1
                    + self.inter_nodes
                    + self.para1_nodes
                    + self.para2_nodes
                )
                / (self.ranvier_nodes - 1)
            )
        return n_comp

    def get_n_segments(self, downsampled):
        """Get number of segments depending on sampling."""
        n_comp = self.get_n_comp(downsampled)
        return (self.n_Ranvier - 1) * n_comp + 1

    def update_axon_morphology(
        self,
        fiber_diam: float,
        axon_length: Optional[float] = None,
        n_Ranvier: Optional[int] = None,
    ) -> dict:
        """Get geometric description of a single axon.

        Parameters
        ----------
        fiber_diam: float
           diameter in micrometers for all fibers in the pathway
        axon_length: float, optional
           axon lengths in mm for all fibers in the pathway
        n_Ranvier: int, optional
           number of nodes of Ranvier per axon.

        Returns
        -------
        dict

        Notes
        -----
        Either axon_length or n_Ranvier needs to be specified.
        """
        self.fiber_diam = fiber_diam
        template = get_axon_parameters_template(fiber_diam)

        # copy from template
        # scaling from mm
        self.ranvier_length = 1e-3 * template["ranvier_length"]
        self.para1_length = 1e-3 * template["para1_length"]
        self.para2_length = 1e-3 * template["para2_length"]
        self.node_step = 1e-3 * template["deltax"]

        # copy without scaling
        self.ranvier_nodes = template["ranvier_nodes"]
        self.inter_nodes = template["inter_nodes"]
        self.para1_nodes = template["para1_nodes"]
        self.para2_nodes = template["para2_nodes"]

        self.n_comp = self.get_n_comp(self.downsampled)

        tmp_inter_length = (
            self.node_step - 2 * self.para1_length - 2 * self.para2_length * 2
        )
        if self.fiber_diam >= 5.7:
            self.inter_length = tmp_inter_length / 6
        else:
            self.inter_length = tmp_inter_length / 3
        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            self.axon_length = axon_length
        else:
            self.n_Ranvier = n_Ranvier

        # additional params for NEURON model, see axon.py
        # TODO make properties
        self.axon_d = template["axon_diameter"]
        self.node_d = template["node_diameter"]
        self.para1_d = template["para1_diameter"]
        self.para2_d = template["para2_diameter"]
        self.lamellas = template["lamellas"]

    # ruff: noqa: C901
    def get_local_compartment_coords(self) -> np.ndarray:
        """Get 1-D coordinates of internodal compartments relative to the node at 0.0.

        Parameters
        ----------
        axon_morphology: dict
            geometric description of a single axon

        Returns
        -------
        Nx1 numpy.ndarray

        """
        loc_coords = np.zeros(self.n_comp - 1, dtype=float)
        loc_pos = 0.0  # just for clarity

        if self.fiber_diam >= 5.7:
            if not self.downsampled:
                # only internodal compartments.
                # The distances will be computed from the node of Ranvier using loc_pos
                for inx_loc in np.arange(1, self.n_comp):
                    inx_loc = int(inx_loc)
                    if inx_loc == 1:
                        loc_pos = (self.ranvier_length + self.para1_length) / 2

                    if inx_loc == 2 or inx_loc == 10:
                        loc_pos = loc_pos + (self.para1_length + self.para2_length) / 2
                    if inx_loc == 3 or inx_loc == 9:
                        loc_pos = loc_pos + (self.para2_length + self.inter_length) / 2
                    if inx_loc in [4, 5, 6, 7, 8]:
                        loc_pos = loc_pos + self.inter_length / 1
                loc_coords[inx_loc - 1] = loc_pos
            else:
                # node -- -- internodal -- -- -- -- internodal -- -- node
                for inx_loc in np.arange(1, self.n_comp):
                    # only internodal compartments.
                    # The distances will be computed from the
                    # node of Ranvier using loc_pos
                    if inx_loc == 1:
                        loc_pos = (
                            (self.ranvier_length + self.inter_length) / 2
                            + self.para1_length
                            + self.para2_length
                        )
                    elif inx_loc == 2:
                        loc_pos = loc_pos + 5 * self.inter_length
                    else:
                        raise RuntimeError("Wrong number of compartments")
                loc_coords[inx_loc - 1] = loc_pos

        elif self.fiber_diam < 5.7:
            if not self.downsampled:
                for inx_loc in np.arange(1, self.n_comp):
                    if inx_loc == 1:
                        loc_pos = (self.ranvier_length + self.para1_length) / 2
                    if inx_loc == 2 or inx_loc == 7:
                        loc_pos = loc_pos + (self.para1_length + self.para2_length) / 2
                    if inx_loc == 3 or inx_loc == 6:
                        loc_pos = loc_pos + (self.para2_length + self.inter_length) / 2
                    if inx_loc == 4 or inx_loc == 5:
                        loc_pos = loc_pos + self.inter_length  # switch to mm from µm
                loc_coords[inx_loc - 1] = loc_pos
            else:
                # mode -- -- -- internodal -- -- -- node
                loc_coords[0] = (
                    0.5 * self.ranvier_length
                    + 1.5 * self.inter_length
                    + self.para1_length
                    + self.para2_length
                )
        return loc_coords


class AxonMorphologyMcNeal1976(AxonMorphology):
    """Axon morphology class for the McNeal1976 model."""

    @property
    def n_comp(self):
        """Only nodes and one internodal per segment."""
        # node -- -- -- internodal -- -- -- node
        return 2

    @property
    def node_step(self):
        """Node step."""
        return self.fiber_diam * 0.2

    @node_step.setter
    def node_step(self, value):
        """Node step."""
        _logger.warning(
            "The node step in the McNeal1976 model is fixed."
            "The value remains unchanged."
        )

    @AxonMorphology.downsampled.setter
    def downsampled(self, value):
        """Downsampled never works for McNeal1976."""
        if value is True:
            raise NotImplementedError("Downsampled McNeal1976 not implemented.")
        self._downsampled = value

    def update_axon_morphology(
        self,
        fiber_diam: float,
        axon_length: Optional[float] = None,
        n_Ranvier: Optional[int] = None,
    ) -> dict:
        """Get geometric description of a single axon.

        Parameters
        ----------
        fiber_diam: float
           diameter in micrometers for all fibers in the pathway
        axon_length: float, optional
           axon lengths in mm for all fibers in the pathway
        n_Ranvier: int, optional
           number of nodes of Ranvier per axon.

        Returns
        -------
        dict

        Notes
        -----
        Either axon_length or n_Ranvier needs to be specified.
        """
        self.fiber_diam = fiber_diam

        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            self.axon_length = axon_length
        else:
            self.n_Ranvier = n_Ranvier

    def get_local_compartment_coords(self):
        """Get 1-D coordinates of internodal compartments relative to the node at 0.0.

        Returns
        -------
        Nx1 numpy.ndarray

        """
        if self.downsampled is True:
            raise NotImplementedError("Downsampled McNeal1976 not implemented.")
        loc_coords = np.zeros(self.n_comp - 1, dtype=float)
        # mode -- -- -- internodal -- -- -- node
        loc_coords[0] = self.node_step * 0.5
        return loc_coords


class AxonModels:
    """Model to represent axons for simulation in OSS-DBS."""

    def __init__(self, stim_dir: str, hemis_idx: int, description_file: str):
        """Model to represent axons for simulation in OSS-DBS
        Fiber trajectories are used to allocate axon models.

        Parameters
        ----------
        stim_dir: str
            full path to the folder where allocated axons are stored
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        description_file: str
            full path to oss-dbs_parameters.mat or a .json file
            that contains the following parameters:
                pathway_mat_file: list
                    full paths to pathways files
                    in lead-dbs format (could be just one)
                axon_diams_all: list
                    diameters in micrometers for all
                    provided fibers, one per pathway
                axon_lengths_all: list
                    axon lengths in mm, one per pathway
                centering_coordinates: list[list]
                    3-D coordinates used to center axons on fibers
                    (e.g. active contacts)
                axon_model: str
                    NEURON model
                combined_h5_file: str
                    full path to the file where axons are stored
                projection_names: list of str, optional
                    Names
                connectome_name: str, optional
                    Connectome name

        Notes
        -----
        oss-dbs_parameters.mat is created via Lead-DBS.
        For .json parameters, see _import_custom_neurons()

        """
        # To find files
        self.stim_dir = stim_dir

        # better safe than sorry because of MATLAB
        hemis_idx = int(hemis_idx)
        if hemis_idx not in [0, 1]:
            raise ValueError("hemis_idx has to be either 0 or 1")

        if hemis_idx == 0:
            self.oss_sim_folder = "OSS_sim_files_rh"
        else:
            self.oss_sim_folder = "OSS_sim_files_lh"

        # defaults, they will be overwritten
        # TODO wrap them in @property ?
        self.pathway_mat_file = None
        self.axon_diams_all = None
        self.axon_lengths_all = None
        self.centering_coordinates = None
        self.projection_names = None
        self.connectome_name = "MyTracts"

        self._read_input(description_file, hemis_idx)

    def _read_input(self, description_file: str, hemis_idx: int):
        """Reads input from input file.

        Notes
        -----
        If a .mat file is provided, it is assumed that
        the information comes from Lead-DBS.
        Otherwise, a custom neuron parser is used.
        """
        _, file_ending = os.path.splitext(description_file)
        if file_ending == ".mat":
            _logger.info("Read Lead-DBS file.")
            self._import_leaddbs_neurons(description_file, hemis_idx)
        elif file_ending == ".json":
            _logger.info("Read custom neuron json file.")
            self._import_custom_neurons(description_file)
        else:
            raise NotImplementedError(
                f"Unsupported input format {file_ending}, "
                "provide either a json or a mat-file."
            )

    @property
    def stim_dir(self):
        """Stimulation directory, where all files are to be found."""
        return self._stim_dir

    @stim_dir.setter
    def stim_dir(self, value: str):
        if not os.path.isdir(value):
            raise ValueError("The provided stimulation directory does not exist.")
        self._stim_dir = value

    @property
    def axon_model(self):
        """Name of the axon model."""
        return self._axon_model

    @axon_model.setter
    def axon_model(self, value):
        valid_models = ["MRG2002", "MRG2002_DS", "McNeal1976"]
        if value not in valid_models:
            ValueError(f"The NEURON model is not valid, use one of {valid_models}.")
        self._axon_model = value

    @property
    def combined_h5_file(self):
        """Name of final HDF5 file."""
        return self._combined_h5_file

    @combined_h5_file.setter
    def combined_h5_file(self, value):
        """Name of final HDF5 file."""
        # must have h5 ending
        _, file_ending = os.path.splitext(value)
        if file_ending != ".h5":
            self._combined_h5_file = value + ".h5"
        else:
            self._combined_h5_file = value

    def _import_leaddbs_neurons(self, description_file, hemis_idx):
        """Import Lead-DBS description for axon models from oss-dbs_parameters.mat.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        description_file: str
            File prepared by Lead-DBS

        """
        # load .mat of different versions (WON'T WORK THIS WAY ATM!)
        try:
            file_inp = h5py.File(description_file, mode="r")
        except ValueError as err:
            raise ValueError(
                "Please, save oss-dbs_parameters using "
                "'save(oss-dbs_parameters_path, 'settings', '-v7.3')'"
            ) from err

        # try to read from .mat
        if "neuronModel" in file_inp["settings"]:
            array_ascii = file_inp["settings"]["neuronModel"][:]
            list_ascii = []
            for i in range(array_ascii.shape[0]):
                list_ascii.append(array_ascii[i][0])
            self.axon_model = "".join(chr(i) for i in list_ascii)
            _logger.debug(f"Use {self.axon_model}")
        else:
            _logger.debug("Use McNeal1976 model by default")
            self.axon_model = "McNeal1976"

        # connectome name within Lead-DBS (e.g. 'Multi-Tract: PetersenLUIC')
        array_ascii = file_inp["settings"]["connectome"][:]
        list_ascii = []
        for i in range(array_ascii.shape[0]):
            list_ascii.append(array_ascii[i][0])
        self.connectome_name = "".join(chr(i) for i in list_ascii)

        # 'Multi-tract' connectomes contain multiple pathways
        # (projections) in separate .mat files
        if "Multi-Tract" in self.connectome_name:
            # this file is pre-filtered connectome assembled in one file in Lead-DBS
            self.pathway_mat_file = [
                os.path.join(
                    self.stim_dir,
                    self.connectome_name.rsplit(" ", 1)[1],
                    f"data{hemis_idx + 1}.mat",
                )
            ]

            self.projection_names = []
            for i in range(len(file_inp["settings"]["connectomeTractNames"][0])):
                ext_string = file_inp[
                    file_inp["settings"]["connectomeTractNames"][0][i]
                ]
                list_ascii = []
                for j in range(ext_string.shape[0]):
                    list_ascii.append(ext_string[j][0])
                projection_name = "".join(chr(i) for i in list_ascii)
                self.projection_names.append(projection_name)
        else:
            self.projection_names = ["default"]
            # this file is pre-filtered connectome in Lead-DBS
            self.pathway_mat_file = [
                os.path.join(
                    self.stim_dir, self.connectome_name, f"data{hemis_idx + 1}.mat"
                )
            ]

        # check which contacts are active to seed axons close to them
        # for StimSets check across all of them
        stimSets = bool(
            file_inp["settings"]["stimSetMode"][0][0]
        )  # if StimSets are used, create a dummy ampl_vector
        if stimSets:
            _logger.info("Use stimSets")
            stim_protocols = np.genfromtxt(
                os.path.join(
                    self.stim_dir,
                    self.oss_sim_folder,
                    f"Current_protocols_{hemis_idx}.csv",
                ),
                dtype=float,
                delimiter=",",
                names=True,
            )

            total_contacts = len(list(stim_protocols[0]))
            total_protocols = stim_protocols.shape[0]

            protocols_array = np.zeros((total_protocols, total_contacts), float)
            ampl_vector = list(stim_protocols[0])  # just initialize

            for j in range(total_protocols):
                protocols_array[j, :] = list(stim_protocols[j])
                for i in range(total_contacts):
                    if not math.isnan(protocols_array[j, i]):
                        ampl_vector[i] = 1.0
        else:
            ampl_vector = list(file_inp["settings"]["Phi_vector"][:, hemis_idx])

        self.centering_coordinates = []
        for i in range(len(ampl_vector)):
            if not (math.isnan(ampl_vector[i])):
                a_ref = file_inp["settings"]["contactLocation"][hemis_idx][0]
                b = file_inp[a_ref]
                self.centering_coordinates.append(b[:, i])

        # hardcoded name for axons pre-filtered by Lead-DBS
        self.combined_h5_file = os.path.join(
            self.stim_dir, self.oss_sim_folder, "Allocated_axons.h5"
        )
        self.output_directory = os.path.dirname(self.combined_h5_file)

        # morphology set in Lead-DBS
        self.axon_lengths_all = list(file_inp["settings"]["axonLength"][:][0][:])
        self.axon_diams_all = list(file_inp["settings"]["fiberDiameter"][:][0][:])

    def _import_custom_neurons(self, description_file):
        """Import custom description for axon models from a .json dictionary.

        Example json input
        custom_dict = {
             'pathway_mat_file': ['SMA_hdp_left.mat',
                                  'SMA_hdp_right.mat',
                                  'gpe2stn_sm_left.mat'],

             # axon diameter and length is the same for all axons within the pathway
             'axon_diams_all': [5.7,5.7,3.0],
             'axon_lengths_all':[20.0,20.0,10.0],

             # in this case, we just have some STN coordinates for left and right in MNI
             'centering_coordinates': [[7.5838, -18.3984, 1.8932],
                                       [-7.5838, -18.3984, 1.8932]],
             'axon_model': 'McNeal1976',
             'combined_h5_file': 'dataset/all_tracts'
         }

        """
        with open(description_file) as fp:
            custom_dict = json.load(fp)

        self.pathway_mat_file = custom_dict["pathway_mat_file"]
        self.axon_diams_all = custom_dict["axon_diams_all"]
        self.axon_lengths_all = custom_dict["axon_lengths_all"]
        self.centering_coordinates = custom_dict["centering_coordinates"]
        self.axon_model = custom_dict["axon_model"]
        self.combined_h5_file = custom_dict["combined_h5_file"]
        self.output_directory = os.path.dirname(self.combined_h5_file)

        if "projection_names" in custom_dict:
            self.projection_names = custom_dict["projection_names"]
        else:
            self.projection_names = [
                pathway_file.rsplit(os.sep, 1)[1][0:-4]
                for pathway_file in self.pathway_mat_file
            ]

        if "connectome_name" in custom_dict:
            self.connectome_name = custom_dict["connectome_name"]

    def _select_axon_morphology_model(self):
        if "MRG2002" in self.axon_model:
            if "MRG2002_DS" == self.axon_model:
                ax_morph_model = AxonMorphologyMRG2002(downsampled=True)
            else:
                ax_morph_model = AxonMorphologyMRG2002()
        else:
            ax_morph_model = AxonMorphologyMcNeal1976()
        return ax_morph_model

    def _get_local_axons_fibers(self, i: int, axon_morphology: AxonMorphology):
        """Get local information."""
        # multiple .mat files (manual input)
        if len(self.pathway_mat_file) > 1:
            return self._deploy_axons_fibers(
                self.pathway_mat_file[i],
                self.projection_names[i],
                axon_morphology,
                False,
            )

        # multiple pathways in one .mat file (Lead-DBS dMRI_MultiTract connectome)
        elif "Multi-Tract" in self.connectome_name:
            return self._deploy_axons_fibers(
                self.pathway_mat_file[0],
                self.projection_names[i],
                axon_morphology,
                True,
            )

        # one .mat file without pathway differentiation (Lead-DBS dMRI connectome)
        else:
            return self._deploy_axons_fibers(
                self.pathway_mat_file[0],
                self.projection_names[i],
                axon_morphology,
                False,
            )

    def convert_fibers_to_axons(self):
        """Seed axons iterating over all pathways."""
        # within a projection (pathway), number of nodes of Ranvier per axon is fixed
        n_Ranvier_per_projection_all = np.zeros(
            shape=len(self.axon_lengths_all), dtype=int
        )
        n_Neurons_all = np.zeros(shape=len(self.axon_lengths_all), dtype=int)
        orig_n_Neurons_all = np.zeros(shape=len(self.axon_lengths_all), dtype=int)

        axon_morphology = self._select_axon_morphology_model()

        # iterate over projections (fibers) and seed axons
        for i in range(len(self.axon_diams_all)):
            # various geometric parameters for a single axon
            axon_morphology.update_axon_morphology(
                self.axon_diams_all[i], self.axon_lengths_all[i]
            )

            (
                n_Ranvier_per_projection,
                n_Neurons,
                orig_n_Neurons,
            ) = self._get_local_axons_fibers(i, axon_morphology)
            n_Ranvier_per_projection_all[i] = n_Ranvier_per_projection
            n_Neurons_all[i] = n_Neurons
            orig_n_Neurons_all[i] = orig_n_Neurons

            _logger.info(
                f"{n_Neurons_all[i]} axons seeded for "
                f"{self.projection_names[i]} with "
                f"{n_Ranvier_per_projection_all[i]}"
                " nodes of Ranvier"
            )

        # only add axon diameters for seeded axons
        self.axon_diams = []
        n_Ranvier_per_projection = []
        n_Neurons = []
        orig_n_Neurons = []
        for i in range(len(self.axon_diams_all)):
            if n_Ranvier_per_projection_all[i] != 0:
                self.axon_diams.append(float(self.axon_diams_all[i]))
                n_Ranvier_per_projection.append(int(n_Ranvier_per_projection_all[i]))
                n_Neurons.append(int(n_Neurons_all[i]))
                orig_n_Neurons.append(int(orig_n_Neurons_all[i]))

        self._save_axon_parameters_in_json(
            n_Ranvier_per_projection, n_Neurons, orig_n_Neurons
        )

    def _save_axon_parameters_in_json(
        self, n_Ranvier_per_projection, n_Neurons, orig_n_Neurons
    ):
        """Save minimally required axon description in a .json file.

        Parameters
        ----------
        n_Ranvier_per_projection: list
            number of nodes of Ranvier for axons of each pathway (one entry per pathway)
        n_Neurons: list
            number of neurons seeded per pathway
        orig_n_Neurons: list
            number of neurons per pathway as defined
            in the connectome (before Kuncel pre-filtering)
        """
        # dictionary to store axon parameters
        axon_dict = {
            "n_Ranvier": n_Ranvier_per_projection,
            "axon_diams": self.axon_diams,
            "Axon_Model_Type": self.axon_model,
            "Name_prepared_neuron_array": self.combined_h5_file,
            "Neuron_model_array_prepared": True,
            "N_seeded_neurons": n_Neurons,
            "N_orig_neurons": orig_n_Neurons,
            "connectome_name": self.connectome_name,
        }

        with open(
            self.output_directory + "/Allocated_axons_parameters.json", "w"
        ) as save_as_dict:
            json.dump(axon_dict, save_as_dict)

    def _deploy_axons_fibers(
        self,
        pathway_file: str,
        projection_name: str,
        axon_morphology: AxonMorphology,
        multiple_projections_per_file: bool = False,
    ):
        """Convert streamlines (fibers) to axons and store in OSS-DBS supported format.

        Parameters
        ----------
        pathway_file: str
            full path to .mat file containing fiber descriptions (Lead-DBS format)
        projection_name: str
            pathway name
        axon_morphology: AxonMorphology
            geometric description of a single axon
        multiple_projections_per_file: bool, optional
            flag if pathway_file contains multiple pathways

        Returns
        -------
        int
            number of nodes of Ranvier for axons in this pathway.
            Returns 0 if failed to see (fiber is too short)
        int
            number of axons seeded for the pathway

        TODO third output missing

        Notes
        -----
        Pathways are stored as separate groups in the specified .h5 file.
        Axons are stored in separate 2-D datasets.
        For Paraview visualization, use axon_array_2D_<projection_name>
        For Lead-DBS visualization, use <projection_name>_axons.mat

        """
        # TODO fallback for non hdf5
        try:
            file = h5py.File(pathway_file, mode="r")
        except ValueError:
            _logger.warning("Fell back to MATLAB file")
            file = scipy.io.loadmat(pathway_file)

        if multiple_projections_per_file is False:
            # fiber_array has 4 columns (x,y,z,fiber_index)
            # rows - all points
            fiber_array = file["fibers"][:]
        else:
            fiber_array = file[projection_name]["fibers"][:]

        if fiber_array.ndim == 1:
            _logger.warning(
                f"{projection_name}: Projection is empty"
                " check settings for fibre diameter and axon length."
                " No nodes were seeded."
            )
            return 0, 0, 0
        else:
            # flip check
            if fiber_array.shape[1] == 4 and fiber_array.shape[0] != 4:
                fiber_array = fiber_array.T
                idx_shape_inx = 0
            else:
                idx_shape_inx = 1

            if multiple_projections_per_file is False:
                if "origNum" in file:
                    orig_N_fibers = int(file["origNum"][0][0])
                else:
                    orig_N_fibers = int(file["idx"][:].shape[idx_shape_inx])
            else:
                if "origNum" in file[projection_name]:
                    orig_N_fibers = int(file[projection_name]["origNum"][0][0])
                else:
                    orig_N_fibers = int(
                        file[projection_name]["idx"][:].shape[idx_shape_inx]
                    )

        # covert fiber table to nibabel streamlines
        streamlines = convert_fibers_to_streamlines(fiber_array)

        # resample streamlines to nodes of Ranvier
        streamlines_resampled, _ = resample_fibers_to_Ranviers(
            streamlines, axon_morphology.node_step, axon_morphology.n_Ranvier
        )

        # truncate streamlines to match selected axon length
        # axons are seeded on the segment closest to active contacts
        # or other ROI, see self.centering_coordinates
        streamlines_axons = place_axons_on_streamlines(
            streamlines_resampled, axon_morphology.n_Ranvier, self.centering_coordinates
        )

        # streamlines_axons already contain the position of Ranvier nodes.
        # Now we get internodal compartments
        # and store all coordinates in a 3D array:
        # compartment index, spatial axis, axon index
        axon_array = np.zeros(
            (axon_morphology.n_segments, 3, len(streamlines_axons)), dtype=float
        )

        # 2-D version for Paraview visualization
        axon_array_2D = np.zeros(
            (axon_morphology.n_segments * len(streamlines_axons), 4), dtype=float
        )

        # save axons as separate datasets within groups that correspond to pathways
        # TODO Why is the h5py file opened in 'a' mode?
        hf = h5py.File(self.combined_h5_file, "a")
        # TODO this part fails if the file already existed.
        g = hf.create_group(projection_name)

        # get local coordinates for internodal compartments
        local_comp_coords = axon_morphology.get_local_compartment_coords()
        glob_ind = 0
        for inx_axn in range(len(streamlines_axons)):
            inx_comp = 0
            for inx in range(axon_morphology.n_Ranvier - 1):
                # compartments are seeded along the internodal vector
                internodal_vector_normalized = normalized(
                    streamlines_axons[inx_axn][inx + 1]
                    - streamlines_axons[inx_axn][inx]
                )

                # positions for nodes of Ranvier are known
                axon_array[inx_comp, :, inx_axn] = streamlines_axons[inx_axn][inx]

                # now place the compartments until the next node
                for loc_comp_inx in range(1, axon_morphology.n_comp):
                    axon_array[inx_comp + loc_comp_inx, :, inx_axn] = (
                        axon_array[inx_comp, :, inx_axn]
                        + local_comp_coords[loc_comp_inx - 1]
                        * internodal_vector_normalized[0][:]
                    )

                inx_comp = inx_comp + axon_morphology.n_comp

            # last node of Ranvier
            axon_array[-1, :, inx_axn] = streamlines_axons[inx_axn][-1]

            axon_array_2D[glob_ind : glob_ind + axon_morphology.n_segments, :3] = (
                axon_array[:, :, inx_axn]
            )
            axon_array_2D[glob_ind : glob_ind + axon_morphology.n_segments, 3] = (
                inx_axn + 1
            )  # because in Matlab they start from 1

            g.create_dataset("axon" + str(inx_axn), data=axon_array[:, :, inx_axn])

            glob_ind = glob_ind + axon_morphology.n_segments

        hf.close()

        return axon_morphology.n_Ranvier, len(streamlines_axons), orig_N_fibers
