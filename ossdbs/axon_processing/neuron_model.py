import fileinput
import json
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from importlib import resources

import h5py
import neuron
import numpy as np
from axon import Axon
from scipy.io import savemat

_logger = logging.getLogger(__name__)


class NeuronSimulator(ABC):
    """Interface to NEURON simulator."""

    # directory, where all NEURON simulations will be conducted
    _neuron_workdir = "neuron_model"
    # by default we don't assume a downsampled model
    _downsampled = False
    # needs to be implemented
    _axon_model = None

    def __init__(self, pathways_dict: dict, output_path: str):
        # create local directory if it does not exist
        if not os.path.isdir(self._neuron_workdir):
            os.mkdir(self._neuron_workdir)

        self._neuron_executable = "nrnivmodl"
        # executable named different on Windows
        if sys.platform == "win32":
            self._neuron_executable = "mknrndll"
        self._pathways_dict = pathways_dict
        self._output_path = output_path  # TODO needs to be full path or relative?

        # create output path if it does not exist
        if not os.path.isdir(output_path):
            os.mkdir(self.output_path)
        # copy NEURON files and compile them
        self.copy_neuron_files()
        self.compile_neuron_files()

        # time-domain solution (not yet loaded)
        self._td_solution = None
        # signal dict will be created when
        # time-domain solution is loaded
        self._signal_dict = None

    @property
    def axon_model(self):
        return self._axon_model

    @property
    def signal_dict(self):
        return self._signal_dict

    @signal_dict.setter
    def signal_dict(self, value: dict):
        self._signal_dict = value

    @property
    def axon_model_type(self):
        return self._pathways_dict["Axon_Model_Type"]

    @property
    def connectome_name(self):
        return self._pathways_dict["connectome_name"]

    def get_pathway_name(self, idx: int):
        return self._pathways_dict["pathway_name"][idx]

    def get_axon_diam(self, idx: int):
        return self._pathways_dict["axon_diams"][idx]

    def get_n_Ranvier(self, idx: int):
        return (self._pathways_dict["n_Ranvier"][idx],)

    def get_N_seeded_neurons(self, idx: int):
        return self._pathways_dict["N_seeded_neurons"][idx]

    def get_N_orig_neurons(self, idx: int):
        return self._pathways_dict["N_orig_neurons"][idx]

    @property
    def downsampled(self) -> bool:
        return self._downsampled

    @downsampled.setter
    def downsampled(self, value: bool) -> None:
        self._downsampled = value

    @property
    def output_path(self):
        """Output path to save results."""
        return self._output_path

    @property
    def resources_path(self):
        """Path were NEURON templates are stored."""
        raise NotImplementedError("Need to add the resource path in source code.")

    @abstractmethod
    def get_axon_morphology(self, axon_diam, axon_length=None, n_Ranvier=None) -> dict:
        """Get geometric description of a single axon.

        Parameters
        ----------
         axon_diam: float
            diameter in micrometers for all fibers in the pathway
         axon_length: float, optional
            axon lengths in mm for all fibers in the pathway. If not specified, provide n_Ranvier
         n_Ranvier: int, optional
            number of nodes of Ranvier per axon. If not specified, provide axon_l
        ength.

        Returns
        -------
        dict

        """

    def copy_neuron_files(self):
        """Copy files from template folder to local folder."""
        path_to_ossdbs = resources.files("ossdbs")
        path_to_neuron_files = path_to_ossdbs.joinpath(
            "axon_processing", "neuron_templates", self.resources_path
        )
        pwd = os.getcwd()
        local_neuron_path = os.path.join(pwd, self._neuron_workdir)
        shutil.copytree(path_to_neuron_files, local_neuron_path)

    @property
    def pathways_dict(self):
        """Information about pathways."""
        return self._pathways_dict

    @abstractmethod
    def modify_hoc_file(self):
        """Update parameters in the hoc file."""
        pass

    @property
    def neuron_executable(self) -> str:
        """Name of the NEURON executable."""
        return self._neuron_executable

    @abstractmethod
    def paste_to_hoc(self):
        """Paste Python parameters into HOC file."""
        pass

    @abstractmethod
    def paste_paraview_vis(self):
        """Convert to paraview file."""

    @abstractmethod
    def compile_neuron_files(self):
        """Compile a NEURON file."""

    def load_solution(self, time_domain_h5_file: str):
        """Load solution from h5 file.

        Parameters
        ----------
        time_domain_h5_file: str
            Name of file holding time-domain solution
        """
        self._td_solution = h5py.File(time_domain_h5_file, "r")

    def create_leaddbs_outputs(self, Axon_Lead_DBS):
        """Export axons with activation state in Lead-DBS supported format.

        Parameters
        ----------
        Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

        """
        mdic = {
            "fibers": Axon_Lead_DBS,
            "ea_fibformat": "1.0",
            "connectome_name": self.connectome_name,
        }  # For Lead-DBS .mat files

        if self.scaling_index is None:
            if self.pathway_name is None:
                savemat(self.output_path + "/Axon_state.mat", mdic)
            else:
                savemat(
                    self.output_path + "/Axon_state_" + self.pathway_name + ".mat", mdic
                )
        else:
            if self.pathway_name is None:
                savemat(
                    self.output_path
                    + "/Axon_state_"
                    + str(self.scaling_index)
                    + ".mat",
                    mdic,
                )
            else:
                savemat(
                    self.output_path
                    + "/Axon_state_"
                    + self.pathway_name
                    + "_"
                    + str(self.scaling_index)
                    + ".mat",
                    mdic,
                )

    def create_paraview_outputs(self, Axon_Lead_DBS):
        """Export axons with activation state in Paraview supported format.

        Parameters
        ----------
        Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

        """
        if self.scaling_index is None:
            if self.pathway_name is None:
                np.savetxt(
                    self.output_path + "/Axon_state.csv",
                    Axon_Lead_DBS,
                    delimiter=",",
                    header="x-pt,y-pt,z-pt,idx,status",
                )
            else:
                np.savetxt(
                    self.output_path + "/Axon_state_" + self.pathway_name + ".csv",
                    Axon_Lead_DBS,
                    delimiter=",",
                    header="x-pt,y-pt,z-pt,idx,status",
                )
        else:
            if self.pathway_name is None:
                np.savetxt(
                    self.output_path
                    + "/Axon_state_"
                    + str(self.scaling_index)
                    + ".csv",
                    Axon_Lead_DBS,
                    delimiter=",",
                    header="x-pt,y-pt,z-pt,idx,status",
                )
            else:
                np.savetxt(
                    self.output_path
                    + "/Axon_state_"
                    + self.pathway_name
                    + "_"
                    + str(self.scaling_index)
                    + ".csv",
                    Axon_Lead_DBS,
                    delimiter=",",
                    header="x-pt,y-pt,z-pt,idx,status",
                )

    def store_axon_statuses(
        self,
        percent_activated,
        percent_damaged,
        percent_csf,
        scaling_index=None,
        pathway_name=None,
    ):
        """Store PAM results.

        Parameters
        ----------
        percent_activated: float, percent of the original(!) number of neurons that are activated for the particular stimulation
        percent_damaged: float, percent of the original(!) number of neurons that are 'damaged' for the particular electrode placement
        percent_csf: float, percent of the original(!) number of neurons that intersect with CSF for the particular brain segmentation

        Note
        ----------
        For activation state of particular neuron see 'fiberActivation*' files as those restore original(!) indices as in the connectome.
        """
        summary_dict = {
            "percent_activated": percent_activated,
            "percent_damaged": percent_damaged,
            "percent_csf": percent_csf,
        }

        if self.scaling_index is None:
            if self.pathway_name is None:
                path_to_save = os.path.join(self.output_path, "Pathway_status.json")
            else:
                summary_dict["pathway_name"] = self.pathway_name
                path_to_save = os.path.join(
                    self.output_path, f"Pathway_status_{pathway_name}.json"
                )
        else:
            summary_dict["scaling_index"] = str(self.scaling_index)
            if self.pathway_name is None:
                path_to_save = os.path.join(
                    self.output_path, f"Pathway_status_{scaling_index}.json"
                )
            else:
                summary_dict["pathway_name"] = self.pathway_name
                path_to_save = os.path.join(
                    self.output_path,
                    f"Pathway_status_{pathway_name}_{scaling_index}.json",
                )
        with open(path_to_save, "w") as save_as_dict:
            json.dump(summary_dict, save_as_dict)

    def process_pathways(self, scaling: float, scaling_index):
        """Go through all pathways and compute the activation.

        Parameters
        ----------
        scaling: float
            scaling factor for the whole solution (different from scaling_vector)
        scaling_index: int
            index of the scaling factor or scaling vector
        """
        if self._td_solution is None:
            raise ValueError("Need to load time-domain solution from H5 file first.")
        pathways = list(self._td_solution.keys())
        pathways.remove("TimeSteps[s]")

        # signal parameters can be extracted from solution
        TimeSteps = np.array(self._td_solution["TimeSteps[s]"])

        self.signal_dict = {
            "time_step": np.round(
                1000.0 * (TimeSteps[1] - TimeSteps[0]), 6
            ),  # convert to ms
            "scaling": scaling,  # from Lead-DBS GUI
            "N_time_steps": TimeSteps.shape[0],
        }

        _logger.info("Going through pathways")
        for pathway_idx, pathway_name in enumerate(pathways):
            pathway_dataset = self._td_solution[pathway_name]
            assert pathway_name == self.get_pathway_name(pathway_idx)
            self.check_pathway_activation(pathway_dataset, pathway_idx)

        # close H5 file
        self._td_solution.close()

    def get_v_ext(self, v_time_sol):
        """Convert el. potential computed by OSS-DBS to extracellular el. potential taking into account scaling parameters.

        Parameters
        ----------
        v_time_sol: NxM numpy.ndarray, (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        v_ext: NxM numpy.ndarray, extracellular el. potential distribution (mV) in space (on the neuron) and time (DBS signal)

        """
        v_ext = np.zeros_like(v_time_sol, float)
        if self.scaling_vector is None:
            v_ext = v_time_sol * 1000.0 * self.scaling  # convert to mV
        else:
            for contact_i in range(len(self.scaling_vector)):
                v_ext = (
                    v_ext
                    + v_time_sol
                    * 1000.0
                    * self.scaling
                    * self.scaling_vector[contact_i]
                )  # convert to mV

        return v_ext

    def get_axon_status_multiprocessing(self, neuron_index, v_time_sol, output):
        """Probe action potential at a neuron (axon).

        Parameters
        ----------
        neuron_index: int
            index of neuron in the pathway starting from 0
        v_time_sol: np.ndarray
            (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)
        output:
            multiprocessing output data structure

        Returns
        -------
        list, neuron index in the pathway and its activation status (1 or 0)

        """
        if self.downsampled:
            v_time_sol = self.upsample_voltage(v_time_sol)
        v_ext = self.get_v_ext(v_time_sol)
        spike = self.run_NEURON(
            self.v_init, v_ext, self.hoc_file, self.extra_initialization
        )
        return output.put([neuron_index, spike])

    def run_NEURON(
        self,
        v_init: float,
        v_ext: np.ndarray,
        hoc_file: str,
        extra_initialization: bool,
    ) -> bool:
        """Load a NEURON model and run simulation.

        Parameters
        ----------
        v_init: float
            Initial potential, resting potential
        v_time_sol: np.ndarray
            NxM extracellular potential distribution in mV (N: neuron, M: time)
        hoc_file: str
            Name of .hoc file to be loaded
        extra_initialization: bool
            Delete and create nodes (needed for MRG2002 model)
        """
        # prints {load_file("...")}
        neuron.h(f'{{load_file("{hoc_file}")}}')
        if extra_initialization:
            neuron.h.deletenodes()
            neuron.h.createnodes()
            neuron.h.dependent_var()
            neuron.h.initialize()

        neuron.h.setupAPWatcher_0()  # 'left' end of axon
        neuron.h.setupAPWatcher_1()  # 'right' end of axon

        neuron.h.dt = self.signal_dict["time_step"]
        neuron.h.tstop = (
            self.signal_dict["N_time_steps"] * self.signal_dict["time_step"]
        )
        # Number of DBS pulses is 1! Hardcoded (TODO)
        neuron.h.n_pulse = 1
        neuron.h.v_init = v_init

        # feed the potential in time for compartment i to the NEURON model
        for i in range(v_ext.shape[0]):
            neuron.h.wf[i] = neuron.h.Vector(v_ext[i, :])

        neuron.h.stimul()
        neuron.h.run()
        # decide if activated
        activated = bool(neuron.h.stoprun)
        return activated

    def check_pathway_activation(self, pathway_dataset, pathway_idx):
        """Parallelized probing of action potentials at all stimulated neurons (axons) described by supplied pathway datagroups.

        Parameters
        ----------
        pathway_dataset: h5 group, contains datasets that describe geometries for all neuron (axons)
        pre_status: Nx1 numpy.ndarray, vector of initial neuron (axon) statuses (0 - available for probing, -1 - intersected with electrode, -2 - inside CSF)

        Note
        ----------
        Creates 'Axon_states_*' .mat and .csv files for visualization in Lead-DBS and Paraview, respectively.
        Also stores summary statistics in 'Pathway_status_*.json'
        """
        # use half of CPUs
        N_proc = mp.cpu_count() / 2

        # get parameters
        N_neurons = self.get_N_seeded_neurons(pathway_idx)
        pathway_name = self.get_pathway_name(pathway_idx)
        axon_diam = self.get_axon_diam(pathway_idx)
        n_Ranvier = self.get_n_Ranvier(pathway_idx)
        orig_N_neurons = self.get_N_orig_neurons(pathway_idx)

        # check actual number of n_segments in case downsampeld
        ax_mh = self.get_axon_morphology(self.axon_model, axon_diam, None, n_Ranvier)
        n_segments_actual = ax_mh["n_segments"]

        # TODO
        self.modify_hoc_file()

        # check pre-status
        pre_status = pathway_dataset["Status"]

        # initialize outputs
        Axon_Lead_DBS = np.zeros((N_neurons * n_segments_actual, 5), float)
        List_of_activated = []
        List_of_not_activated = []
        Activated_models = 0

        neuron_index = 0
        while neuron_index < N_neurons:
            # run parallel processing
            proc = []
            j_proc = 0  # counter for processes
            output = mp.Queue()
            while j_proc < N_proc and neuron_index < N_neurons:
                # get neuron geometry and field solution
                neuron = pathway_dataset["axon" + str(neuron_index)]
                Axon_Lead_DBS[
                    neuron_index
                    * self.n_segments_actual : (neuron_index + 1)
                    * self.n_segments_actual,
                    :3,
                ] = np.array(neuron["Points[mm]"])

                # add index
                Axon_Lead_DBS[
                    neuron_index
                    * self.n_segments_actual : (neuron_index + 1)
                    * self.n_segments_actual,
                    3,
                ] = (
                    neuron_index + 1
                )  # because Lead-DBS numbering starts from 1

                # check which neurons were flagged with CSF and electrode intersection, skip probing of those
                if pre_status[neuron_index] != 0:
                    Axon_Lead_DBS[
                        neuron_index
                        * self.n_segments_actual : (neuron_index + 1)
                        * self.n_segments_actual,
                        4,
                    ] = pre_status[neuron_index]
                    neuron_index += 1
                    continue

                neuron_time_sol = np.array(neuron["Potential[V]"])

                processes = mp.Process(
                    target=self.get_axon_status,
                    args=(neuron_index, neuron_time_sol, output),
                )
                proc.append(processes)

                j_proc += 1
                neuron_index += 1

            for p in proc:
                p.start()
            for p in proc:
                p.join()

            # check the status of batch processed neurons
            neurons_idxs_stat = [output.get() for p in proc]
            # n_idx_stat is a list[neuron index, status (1 or 0)]
            for n_idx_stat in neurons_idxs_stat:
                if n_idx_stat[1] == 1:
                    Activated_models += 1
                    List_of_activated.append(n_idx_stat[0])
                else:
                    List_of_not_activated.append(n_idx_stat[0])

        # iterate over all neurons initially placed by OSS-DBS and assign their statuses
        for neuron_index in range(N_neurons):
            if neuron_index in List_of_activated:
                Axon_Lead_DBS[
                    neuron_index
                    * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    4,
                ] = 1
            elif neuron_index in List_of_not_activated:
                Axon_Lead_DBS[
                    neuron_index
                    * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    4,
                ] = 0
            else:
                # the status was already assigned
                continue

        self.create_leaddbs_outputs(Axon_Lead_DBS)
        self.create_paraview_outputs(Axon_Lead_DBS)

        percent_activated = np.round(
            100.0 * Activated_models / float(orig_N_neurons), 2
        )
        percent_damaged = np.round(
            100.0 * np.sum(np.isclose(pre_status, -1.0)) / float(orig_N_neurons), 2
        )
        percent_csf = np.round(
            100.0 * np.sum(np.isclose(pre_status, -2.0)) / float(orig_N_neurons), 2
        )

        _logger.info(f"\n\nPathway {pathway_name} :")
        _logger.info(f"Activated neurons: {percent_activated}%")
        _logger.info(f"Neurons damaged: {percent_damaged}%")
        _logger.info(f"Neurons in CSF {percent_csf}%")

        self.store_axon_statuses(percent_activated, percent_damaged, percent_csf)


class MRG2002(NeuronSimulator):
    _axon_model = "MRG2002"

    @property
    def resources_path(self):
        return "MRG2002"

    def get_axon_morphology(self, axon_diam, axon_length=None, n_Ranvier=None) -> dict:
        """Get geometric description of a single axon.

        Parameters
        ----------
         axon_diam: float
            diameter in micrometers for all fibers in the pathway
         axon_length: float, optional
            axon lengths in mm for all fibers in the pathway. If not specified, provide n_Ranvier
         n_Ranvier: int, optional
            number of nodes of Ranvier per axon. If not specified, provide axon_l
        ength.

        Returns
        -------
        dict


        TODO rewrite

        """
        axon_morphology = {"axon_model": self.axon_model, "axon_diam": axon_diam}
        param_ax = {"centered": True, "diameter": axon_diam}
        a = Axon(param_ax)
        nr = Axon.get_axonparams(a)

        axon_morphology["ranvier_length"] = 1e-3 * nr["ranvier_length"]
        axon_morphology["para1_length"] = 1e-3 * nr["para1_length"]
        axon_morphology["para2_length"] = 1e-3 * nr["para2_length"]
        axon_morphology["node_step"] = 1e-3 * nr["deltax"]
        if self.downsampled:
            if axon_diam >= 5.7:
                # node -- -- internodal -- -- -- -- internodal -- -- node
                axon_morphology["n_comp"] = 3
                axon_morphology["inter_length"] = (
                    axon_morphology["node_step"]
                    - axon_morphology["para1_length"] * 2
                    - axon_morphology["para2_length"] * 2
                ) / 6
            else:
                # node -- -- -- internodal -- -- -- node
                axon_morphology["n_comp"] = 2
                axon_morphology["inter_length"] = (
                    axon_morphology["node_step"]
                    - axon_morphology["para1_length"] * 2
                    - axon_morphology["para2_length"] * 2
                ) / 3
        else:
            axon_morphology["n_comp"] = int(
                (
                    (nr["ranvier_nodes"] - 1)
                    + nr["inter_nodes"]
                    + nr["para1_nodes"]
                    + nr["para2_nodes"]
                )
                / (nr["ranvier_nodes"] - 1)
            )

            if axon_diam >= 5.7:
                axon_morphology["inter_length"] = (
                    axon_morphology["node_step"]
                    - axon_morphology["para1_length"] * 2
                    - axon_morphology["para2_length"] * 2
                ) / 6
            else:
                axon_morphology["inter_length"] = (
                    axon_morphology["node_step"]
                    - axon_morphology["para1_length"] * 2
                    - axon_morphology["para2_length"] * 2
                ) / 3
        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            axon_morphology["axon_length"] = axon_length
            axon_morphology["n_Ranvier"] = int(
                axon_length / axon_morphology["node_step"]
            )
        else:
            axon_morphology["n_Ranvier"] = n_Ranvier
            axon_morphology["axon_length"] = n_Ranvier * axon_morphology["node_step"]

        # always odd number of nodes of Ranvier!
        if axon_morphology["n_Ranvier"] % 2 == 0:
            axon_morphology["n_Ranvier"] -= 1
            axon_morphology["axon_length"] = (
                axon_morphology["n_Ranvier"] * axon_morphology["node_step"]
            )

        axon_morphology["n_para1"] = (
            nr["para1_nodes"] * (axon_morphology["n_Ranvier"] - 1) / (21 - 1)
        )
        axon_morphology["n_para2"] = (
            nr["para2_nodes"] * (axon_morphology["n_Ranvier"] - 1) / (21 - 1)
        )

        axon_morphology["n_segments"] = int(
            (axon_morphology["n_Ranvier"] - 1) * axon_morphology["n_comp"] + 1
        )  # overall number of points on Axon incl. internodal

        # additional params for NEURON model, see axon.py
        (
            axon_morphology["axon_d"],
            axon_morphology["node_d"],
            axon_morphology["para1_d"],
            axon_morphology["para2_d"],
            axon_morphology["lamellas"],
        ) = (
            nr["axon_diameter"],
            nr["node_diameter"],
            nr["para1_diameter"],
            nr["para2_diameter"],
            nr["lamellas"],
        )
        return axon_morphology

    def paste_to_hoc(
        self,
        axonnodes,
        paranodes1,
        paranodes2,
        axoninter,
        axontotal,
        v_init,
        fiberD,
        paralength1,
        paralength2,
        nodelength,
        nodeD,
        axonD,
        paraD1,
        paraD2,
        deltax,
        nl,
        steps_per_ms,
    ):
        axonnodes_line = "axonnodes="
        axonnodes_input = f"axonnodes={axonnodes}\n"

        paranodes1_line = "paranodes1="
        paranodes1_input = f"paranodes1={paranodes1}\n"

        paranodes2_line = "paranodes2="
        paranodes2_input = f"paranodes2={paranodes2}\n"

        axoninter_line = "axoninter="
        axoninter_input = f"axoninter={axoninter}\n"

        axontotal_line = "axontotal="
        axontotal_input = f"axontotal={axontotal}\n"

        nv_init_line = "v_init="
        nv_init_input = f"v_init={v_init}\n"  # normally, -80mv

        fiberD_line = "fiberD="
        fiberD_input = f"fiberD={fiberD}\n"  # fiber diameter

        paralength1_line = "paralength1="
        paralength1_input = f"paralength1={paralength1}\n"

        paralength2_line = "paralength2="
        paralength2_input = f"paralength2={paralength2}\n"

        nodelength_line = "nodelength="
        nodelength_input = f"nodelength={nodelength}\n"

        nodeD_line = "nodeD="
        nodeD_input = f"nodeD={nodeD}\n"

        axonD_line = "axonD="
        axonD_input = f"axonD={axonD}\n"

        paraD1_line = "paraD1="
        paraD1_input = f"paraD1={paraD1}\n"

        paraD2_line = "paraD2="
        paraD2_input = f"paraD2={paraD2}\n"

        deltax_line = "deltax="
        deltax_input = f"deltax={deltax}\n"

        nl_line = "nl="
        nl_input = f"nl={nl}\n"

        steps_per_ms_line = "steps_per_ms="
        steps_per_ms_input = f"steps_per_ms={steps_per_ms}\n"

        x = fileinput.input(files="axon4pyfull.hoc", inplace=1)
        for line in x:
            if line.startswith(axonnodes_line):
                line = axonnodes_input
            if line.startswith(paranodes1_line):
                line = paranodes1_input
            if line.startswith(paranodes2_line):
                line = paranodes2_input
            if line.startswith(axoninter_line):
                line = axoninter_input
            if line.startswith(axontotal_line):
                line = axontotal_input
            if line.startswith(nv_init_line):
                line = nv_init_input
            if line.startswith(fiberD_line):
                line = fiberD_input
            if line.startswith(paralength1_line):
                line = paralength1_input
            if line.startswith(paralength2_line):
                line = paralength2_input
            if line.startswith(nodelength_line):
                line = nodelength_input

            if line.startswith(nodeD_line):
                line = nodeD_input
            if line.startswith(axonD_line):
                line = axonD_input
            if line.startswith(paraD1_line):
                line = paraD1_input
            if line.startswith(paraD2_line):
                line = paraD2_input
            if line.startswith(deltax_line):
                line = deltax_input
            if line.startswith(nl_line):
                line = nl_input
            if line.startswith(steps_per_ms_line):
                line = steps_per_ms_input
            print(line, end="")
        x.close()

        return True

    def compile_neuron_files(self):
        subprocess.run(
            "nocmodl axnode.mod",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            cwd=self._neuron_workdir,
        )  # might not work with remote hard drives
        # TODO test if works
        subprocess.run(
            self.neuron_executable,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            cwd=self._neuron_workdir,
        )

    def paste_paraview_vis(self, Points_on_model, N_comp_in_between):
        NPoints_line = "Points_on_model"
        NPoints_input = f"Points_on_model={Points_on_model}\n"  # NEURON uses ms
        N_comp_in_between_line = "N_comp_in_between"
        N_comp_in_between_input = (
            f"N_comp_in_between={N_comp_in_between}\n"  # NEURON uses ms
        )

        fl = fileinput.input(
            files="Visualization_files/Paraview_vis_axon.py", inplace=1
        )
        for line in fl:
            if line.startswith(NPoints_line):
                line = NPoints_input
            if line.startswith(N_comp_in_between_line):
                line = N_comp_in_between_input
            print(line, end="")
        fl.close()

        return True

    def modify_hoc_file(self, nRanvier, stepsPerMs, axon_morphology, axonDiam):
        if self.axonDiam >= 5.7:
            axoninter = (nRanvier - 1) * 6
        else:
            axoninter = (nRanvier - 1) * 3

        v_init = -80.0

        self.paste_to_hoc(
            nRanvier,
            axon_morphology["n_para1"],
            axon_morphology["n_para2"],
            axoninter,
            axon_morphology["n_segments"],
            v_init,
            axonDiam,
            axon_morphology["para1_length"] * 1e3,
            axon_morphology["para2_length"] * 1e3,
            axon_morphology["ranvier_length"] * 1e3,
            axon_morphology["node_d"],
            axon_morphology["axon_d"],
            axon_morphology["para1_d"],
            axon_morphology["para2_d"],
            axon_morphology["node_step"] * 1e3,
            axon_morphology["lamellas"],
            stepsPerMs,
        )

    def upsample_voltage(self, v_time_sol):
        """Upsample potential distribution for downsampled neurons by interpolating.

        Parameters
        ----------
        v_time_sol: LxM numpy.ndarray, (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        v_time_sol_full: NxM numpy.ndarray

        Notes
        -----
        TODO estimate ratios directly from the morphology
        """
        # let's interpolate voltage between node - center_l - center_r - node
        # assume 11 segments
        # n_segments_ds = ((n_segments_full - 1) / 11) * 3 +1

        v_time_sol_full = np.zeros(
            (self.axon_morphology["n_segments"], v_time_sol.shape[1]), float
        )

        if self.axonDiam >= 5.7:
            # fill out nodes first
            for k in np.arange(0, self.axon_morphology["n_segments"], 11):
                z = int(k / 11) * 3
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # now two segments in between
            for k in np.arange(3, self.axon_morphology["n_segments"], 11):
                z = int(k / 11) * 3 + 1
                v_time_sol_full[k, :] = v_time_sol[z, :]

            for k in np.arange(8, self.axon_morphology["n_segments"], 11):
                z = int(k / 11) * 3 + 2
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # node -- -- intern -- -- -- -- intern -- -- node  ->  node-para1-para2-intern-intern-intern-intern-intern-intern-para2-para1-node
            internodal_length = (
                self.axon_morphology["node_step"]
                - self.axon_morphology["ranvier_length"]
                - (
                    self.axon_morphology["para2_length"]
                    + self.axon_morphology["para1_length"]
                )
                * 2.0
            ) / 6.0
            dist_node_internode = (
                (self.axon_morphology["ranvier_length"] + internodal_length) * 0.5
                + self.axon_morphology["para1_length"]
                + self.axon_morphology["para2_length"]
            )
            ratio_1 = (
                (
                    self.axon_morphology["ranvier_length"]
                    + self.axon_morphology["para1_length"]
                )
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (
                    self.axon_morphology["para1_length"]
                    + self.axon_morphology["para2_length"]
                )
                * 0.5
                / dist_node_internode
            )
            ratio_3 = internodal_length / dist_node_internode
            ratio_4 = ratio_3 + internodal_length / dist_node_internode

            # now interpolate to the rest
            list_interp = [
                [1, 2],
                [4, 5, 6, 7],
                [9, 10],
            ]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, self.axon_morphology["n_segments"] - 1, 11):
                    if interv == 0:
                        v_time_sol_full[j + 1, :] = (1 - ratio_1) * v_time_sol_full[
                            j, :
                        ] + (ratio_1) * v_time_sol_full[j + 3, :]
                        v_time_sol_full[j + 2, :] = (1 - ratio_2) * v_time_sol_full[
                            j, :
                        ] + (ratio_2) * v_time_sol_full[j + 3, :]
                    elif interv == 1:
                        v_time_sol_full[j + 4, :] = (1 - ratio_3) * v_time_sol_full[
                            j + 3, :
                        ] + ratio_3 * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 5, :] = (1 - ratio_4) * v_time_sol_full[
                            j + 3, :
                        ] + ratio_4 * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 6, :] = (
                            ratio_4 * v_time_sol_full[j + 3, :]
                            + (1 - ratio_4) * v_time_sol_full[j + 8, :]
                        )
                        v_time_sol_full[j + 7, :] = (
                            ratio_3 * v_time_sol_full[j + 3, :]
                            + (1 - ratio_3) * v_time_sol_full[j + 8, :]
                        )
                    else:
                        v_time_sol_full[j + 9, :] = (
                            ratio_2 * v_time_sol_full[j + 8, :]
                            + (1 - ratio_2) * v_time_sol_full[j + 11, :]
                        )
                        v_time_sol_full[j + 10, :] = (
                            ratio_1 * v_time_sol_full[j + 8, :]
                            + (1 - ratio_1) * v_time_sol_full[j + 11, :]
                        )
        else:
            # let's interpolate voltage between node - center - node
            # assume 8 segments
            # fill out nodes first
            for k in np.arange(0, self.axon_morphology["n_segments"], 8):
                z = int(k / 8) * 2
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # now the center between nodes
            for k in np.arange(4, self.axon_morphology["n_segments"], 8):
                z = int(k / 8) * 2 + 1
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # node -- -- -- internodal -- -- -- node  ->  node-para1-para2-intern-intern-intern-para2-para1-node
            dist_node_internode = self.axon_morphology["node_step"] / 2.0
            internodal_length = (
                self.axon_morphology["node_step"]
                - self.axon_morphology["ranvier_length"]
                - (
                    self.axon_morphology["para2_length"]
                    + self.axon_morphology["para1_length"]
                )
                * 2
            ) / 3.0
            ratio_1 = (
                (
                    self.axon_morphology["ranvier_length"]
                    + self.axon_morphology["para1_length"]
                )
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (
                    self.axon_morphology["para1_length"]
                    + self.axon_morphology["para2_length"]
                )
                * 0.5
                / dist_node_internode
            )
            ratio_3 = (
                ratio_2
                + (self.axon_morphology["para2_length"] + internodal_length)
                * 0.5
                / dist_node_internode
            )

            # now interpolate to the rest
            list_interp = [
                [1, 2, 3],
                [5, 6, 7],
            ]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, self.axon_morphology["n_segments"] - 1, 8):
                    if interv == 0:  # ratios based on intercompartment distances
                        v_time_sol_full[j + 1, :] = (1 - ratio_1) * v_time_sol_full[
                            j, :
                        ] + (ratio_1) * v_time_sol_full[j + 4, :]
                        v_time_sol_full[j + 2, :] = (1 - ratio_2) * v_time_sol_full[
                            j, :
                        ] + (ratio_2) * v_time_sol_full[j + 4, :]
                        v_time_sol_full[j + 3, :] = (1 - ratio_3) * v_time_sol_full[
                            j, :
                        ] + (ratio_3) * v_time_sol_full[j + 4, :]
                    else:
                        v_time_sol_full[j + 5, :] = (ratio_3) * v_time_sol_full[
                            j + 4, :
                        ] + (1 - ratio_3) * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 6, :] = (ratio_2) * v_time_sol_full[
                            j + 4, :
                        ] + (1 - ratio_2) * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 7, :] = (ratio_1) * v_time_sol_full[
                            j + 4, :
                        ] + (1 - ratio_1) * v_time_sol_full[j + 8, :]

        return v_time_sol_full


class McNeal1976(NeuronSimulator):
    _axon_model = "McNeal1976"

    @property
    def resources_path(self):
        return "McNeal1976"

    def compile_neuron_files(self):
        subprocess.run(
            self.neuron_executable,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            cwd=self._neuron_workdir,
        )

    def get_axon_morphology(
        self, axon_model, axon_diam, axon_length=None, n_Ranvier=None
    ) -> dict:
        """Get geometric description of a single axon.

        Parameters
        ----------
         axon_model: str
            NEURON model ('MRG2002', 'MRG2002_DS' (downsampled), 'McNeal1976' (classic McNeal's))
         axon_diam: float
            diameter in micrometers for all fibers in the pathway
         axon_length: float, optional
            axon lengths in mm for all fibers in the pathway. If not specified, provide n_Ranvier
         n_Ranvier: int, optional
            number of nodes of Ranvier per axon. If not specified, provide axon_l
        ength.

        Returns
        -------
        dict


        TODO rewrite

        """
        axon_morphology = {"axon_model": axon_model, "axon_diam": axon_diam}

        # node -- -- -- internodal -- -- -- node
        axon_morphology["n_comp"] = 2  # only nodes and one internodal per segment
        axon_morphology["node_step"] = axon_diam * 0.2  # from 1 to 2 mm
        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            axon_morphology["axon_length"] = axon_length
            axon_morphology["n_Ranvier"] = int(
                axon_length / axon_morphology["node_step"]
            )
        else:
            axon_morphology["n_Ranvier"] = n_Ranvier
            axon_morphology["axon_length"] = n_Ranvier * axon_morphology["node_step"]

        # always odd number of nodes of Ranvier!
        if axon_morphology["n_Ranvier"] % 2 == 0:
            axon_morphology["n_Ranvier"] -= 1
            axon_morphology["axon_length"] = (
                axon_morphology["n_Ranvier"] * axon_morphology["node_step"]
            )

        axon_morphology["n_segments"] = int(
            (axon_morphology["n_Ranvier"] - 1) * axon_morphology["n_comp"] + 1
        )
        return axon_morphology

    def paste_to_hoc(self, axonnodes, axoninter, axontotal, v_init, steps_per_ms):
        NNODES_line = "NNODES ="
        NNODES_input = f"NNODES = {axonnodes}\n"

        axonnodes_line = "axonnodes="
        axonnodes_input = f"axonnodes={axonnodes}\n"

        nv_init_line = "v_init="
        nv_init_input = f"v_init={v_init}\n"

        steps_per_ms_line = "steps_per_ms="
        steps_per_ms_input = f"steps_per_ms={steps_per_ms}\n"

        x = fileinput.input(files="init_B5_extracellular.hoc", inplace=1)
        for line in x:
            if line.startswith(axonnodes_line):
                line = axonnodes_input
            if line.startswith(nv_init_line):
                line = nv_init_input
            if line.startswith(steps_per_ms_line):
                line = steps_per_ms_input
            print(line, end="")
        x.close()

        x = fileinput.input(files="axon5.hoc", inplace=1)
        for line in x:
            if line.startswith(NNODES_line):
                line = NNODES_input
            print(line, end="")
        x.close()

        return True

    def modify_hoc_file(self, nRanvier, axon_morphology, stepsPerMs):
        n_internodal = axon_morphology["n_segments"] - nRanvier
        v_init = -70.0
        self.paste_to_hoc(
            nRanvier,
            n_internodal,
            axon_morphology["n_segments"],
            v_init,
            stepsPerMs,
        )

    def paste_paraview_vis(self, Points_on_model, N_comp_in_between):
        # strings to be replaced
        NPoints_line = "Points_on_model"
        NPoints_input = f"Points_on_model={Points_on_model}\n"  # NEURON uses ms
        N_comp_in_between_line = "N_comp_in_between"
        N_comp_in_between_input = (
            f"N_comp_in_between={N_comp_in_between}\n"  # NEURON uses ms
        )

        # TODO why are these pathes hardcorded?
        fl = fileinput.input(
            files="Visualization_files/Paraview_vis_axon.py", inplace=1
        )
        for line in fl:
            if line.startswith(NPoints_line):
                line = NPoints_input
            if line.startswith(N_comp_in_between_line):
                line = N_comp_in_between_input
            print(line, end="")
        fl.close()

        return True
