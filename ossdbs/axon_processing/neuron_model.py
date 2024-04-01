import fileinput
import logging
import multiprocessing as mp
import os
import platform
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from typing import Optional

import h5py
import numpy as np

subversion = platform.python_version_tuple()
if int(subversion[1]) == 8:
    from importlib_resources import files
else:
    from importlib.resources import files

from .axon import Axon
from .utilities import (
    create_leaddbs_outputs,
    create_paraview_outputs,
    store_axon_statuses,
)

_logger = logging.getLogger(__name__)


class NeuronSimulator(ABC):
    """Interface to NEURON simulator."""

    # directory, where all NEURON simulations will be conducted
    _neuron_workdir = os.getcwd()  # "neuron_model"
    # by default we don't assume a downsampled model
    _downsampled = False
    # needs to be implemented
    _axon_model = None
    # Extra steps in NEURON
    _extra_initialization = False
    # Initial voltage
    _v_init = -70.0
    # hoc file for execution
    _hoc_file = None

    def __init__(self, pathways_dict: dict, output_path: str):
        # create local directory if it does not exist
        if not os.path.isdir(self._neuron_workdir):
            os.mkdir(self._neuron_workdir)

        self._neuron_executable = "nrnivmodl"
        # executable named different on Windows
        if sys.platform == "win32":
            self._neuron_executable = "mknrndll"
        self._pathways_dict = pathways_dict
        self._output_path = output_path

        # create output path if it does not exist
        if not os.path.isdir(output_path):
            _logger.info(f"Created directory {output_path}")
            os.mkdir(self.output_path)
        # copy NEURON files and compile them
        _logger.info("Copy and compile NEURON files")
        self.copy_neuron_files()
        self.compile_neuron_files()

        # time-domain solution (not yet loaded)
        self._td_solution = None
        # signal dict will be created when
        # time-domain solution is loaded
        self._signal_dict = None

    @property
    def hoc_file(self):
        if self._hoc_file is None:
            raise NotImplementedError("Need to implement default hoc file!")
        return self._hoc_file

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

    def get_axon_diam(self, idx: int):
        return self._pathways_dict["axon_diams"][idx]

    def get_n_Ranvier(self, idx: int):
        return self._pathways_dict["n_Ranvier"][idx]

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
    def get_axon_morphology(
        self, axon_diam, axon_length=None, n_Ranvier=None, downsampled=False
    ) -> dict:
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
        path_to_ossdbs = files("ossdbs")
        path_to_neuron_files = path_to_ossdbs.joinpath(
            "axon_processing", "neuron_templates", self.resources_path
        )
        pwd = os.getcwd()
        local_neuron_path = os.path.join(pwd, self._neuron_workdir)
        # Directory can already exist
        shutil.copytree(path_to_neuron_files, local_neuron_path, dirs_exist_ok=True)

    @property
    def pathways_dict(self):
        """Information about pathways."""
        return self._pathways_dict

    @abstractmethod
    def modify_hoc_file(self, nRanvier, stepsPerMs, axon_morphology):
        """Update parameters in the hoc file."""
        pass

    @property
    def neuron_executable(self) -> str:
        """Name of the NEURON executable."""
        return self._neuron_executable

    @abstractmethod
    def paste_to_hoc(self, parameters_dict: dict):
        """Paste Python parameters into HOC file."""
        pass

    # TODO needed?
    # @abstractmethod
    # def paste_paraview_vis(self):
    #     """Convert to paraview file."""

    def compile_neuron_files(self):
        """Compile a NEURON file."""
        _logger.info("Compile NEURON executable")
        subprocess.run(
            self.neuron_executable,
            # TODO pipe output
            # stdout=subprocess.STDOUT,
            # stdout=subprocess.DEVNULL,
            # stderr=subprocess.STDOUT,
            cwd=os.path.abspath(self._neuron_workdir),
        )

    def load_solution(self, time_domain_h5_file: str):
        """Load solution from h5 file.

        Parameters
        ----------
        time_domain_h5_file: str
            Name of file holding time-domain solution
        """
        self._td_solution = h5py.File(time_domain_h5_file, "r")

    def process_pathways(
        self, scaling: float = 1.0, scaling_index: Optional[int] = None
    ):
        """Go through all pathways and compute the activation.

        Parameters
        ----------
        scaling: float
            scaling factor for the whole solution
        scaling_index: int
            index of the scaling factor or scaling vector


        Notes
        -----
        TODO Do scaling outside this routine?
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
            self.check_pathway_activation(
                pathway_dataset, pathway_idx, pathway_name, scaling_index
            )

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
        scaling = self.signal_dict["scaling"]
        v_ext = np.zeros_like(v_time_sol, float)
        v_ext = v_time_sol * 1000.0 * scaling  # convert to mV

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
        v_ext = self.get_v_ext(v_time_sol)
        spike = self.run_NEURON(v_ext, self._extra_initialization)
        return output.put([neuron_index, spike])

    def run_NEURON(
        self,
        v_ext: np.ndarray,
        extra_initialization: bool = False,
    ) -> bool:
        """Load a NEURON model and run simulation.

        Parameters
        ----------
        v_ext: np.ndarray
            NxM extracellular potential distribution in mV (N: neuron, M: time)
        extra_initialization: bool
            Delete and create nodes (needed for MRG2002 model)
        """
        _logger.debug(f"Load file: {self.hoc_file}")
        # Do not move the import! Causes trouble
        import neuron

        neuron.h(f'{{load_file("{self.hoc_file}")}}')
        if extra_initialization:
            neuron.h.deletenodes()
            neuron.h.createnodes()
            neuron.h.dependent_var()
            neuron.h.initialize()

        neuron.h.setupAPWatcher_0()  # 'left' end of axon
        neuron.h.setupAPWatcher_1()  # 'right' end of axon

        time_step = self.signal_dict["time_step"]
        run_time = self.signal_dict["N_time_steps"] * self.signal_dict["time_step"]
        _logger.debug(f"Set time step to {time_step} ms")
        _logger.debug(f"Number of time steps is {self.signal_dict['N_time_steps']}")
        _logger.debug(f"Will run neuron for {run_time} ms")
        neuron.h.dt = time_step
        neuron.h.tstop = run_time
        # Number of DBS pulses is 1! Hardcoded (TODO)
        neuron.h.n_pulse = 1
        # Initial potential is hardcoded (TODO?)
        neuron.h.v_init = self._v_init

        # feed the potential in time for compartment i to the NEURON model
        for i in range(v_ext.shape[0]):
            neuron.h.wf[i] = neuron.h.Vector(v_ext[i, :])

        _logger.debug("Stimulation")
        neuron.h.stimul()
        neuron.h.run()
        # decide if activated
        activated = bool(neuron.h.stoprun)
        return activated

    def check_pathway_activation(
        self, pathway_dataset, pathway_idx, pathway_name=None, scaling_index=None
    ):
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
        axon_diam = self.get_axon_diam(pathway_idx)
        n_Ranvier = self.get_n_Ranvier(pathway_idx)
        orig_N_neurons = self.get_N_orig_neurons(pathway_idx)

        # TODO refactor, very confusing to have 2 things being almost the same
        axon_morphology = self.get_axon_morphology(axon_diam, n_Ranvier=n_Ranvier)
        # check actual number of n_segments in case downsampled
        ax_mh = self.get_axon_morphology(
            axon_diam, n_Ranvier=n_Ranvier, downsampled=self.downsampled
        )
        n_segments_actual = ax_mh["n_segments"]
        _logger.debug(f"n_segments_actual: {n_segments_actual}")
        stepsPerMs = int(1.0 / self.signal_dict["time_step"])
        # edit hoc file locally to change parameters
        self.modify_hoc_file(n_Ranvier, stepsPerMs, axon_morphology)

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
                    * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    :3,
                ] = np.array(neuron["Points[mm]"])

                # add index
                Axon_Lead_DBS[
                    neuron_index
                    * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    3,
                ] = (
                    neuron_index + 1
                )  # because Lead-DBS numbering starts from 1

                # check which neurons were flagged with CSF and electrode intersection, skip probing of those
                if pre_status[neuron_index] != 0:
                    Axon_Lead_DBS[
                        neuron_index
                        * n_segments_actual : (neuron_index + 1)
                        * n_segments_actual,
                        4,
                    ] = pre_status[neuron_index]
                    neuron_index += 1
                    continue

                neuron_time_sol = np.array(neuron["Potential[V]"])
                # upsample
                if self.downsampled:
                    _logger.debug(f"Before upsampling: {neuron_time_sol.shape}")
                    neuron_time_sol = self.upsample_voltage(
                        neuron_time_sol, axon_diam, axon_morphology
                    )
                    _logger.debug(f"After upsampling: {neuron_time_sol.shape}")
                processes = mp.Process(
                    target=self.get_axon_status_multiprocessing,
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

        create_leaddbs_outputs(
            self.output_path,
            Axon_Lead_DBS,
            self.connectome_name,
            scaling_index=scaling_index,
            pathway_name=pathway_name,
        )
        create_paraview_outputs(
            self.output_path,
            Axon_Lead_DBS,
            scaling_index=scaling_index,
            pathway_name=pathway_name,
        )

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

        store_axon_statuses(
            self.output_path,
            percent_activated,
            percent_damaged,
            percent_csf,
            scaling_index=scaling_index,
            pathway_name=pathway_name,
        )


class MRG2002(NeuronSimulator):
    _axon_model = "MRG2002"
    # needs extra compilation steps
    _extra_initialization = True
    # different initial resting potential
    _v_init = -80.0
    _hoc_file = "axon4pyfull.hoc"

    @property
    def resources_path(self):
        return "MRG2002"

    def get_axon_morphology(
        self, axon_diam, axon_length=None, n_Ranvier=None, downsampled=False
    ) -> dict:
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
        if downsampled:
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

    def paste_to_hoc(self, parameters_dict: dict) -> None:
        info_to_update = [
            "axonnodes",
            "paranodes1",
            "paranodes2",
            "axoninter",
            "axontotal",
            "v_init",
            "fiberD",
            "paralength1",
            "paralength2",
            "nodelength",
            "nodeD",
            "axonD",
            "paraD1",
            "paraD2",
            "deltax",
            "nl",
            "steps_per_ms",
        ]

        # check that all parameters are in dict
        if not set(info_to_update) == set(parameters_dict.keys()):
            raise ValueError("Need to provide all parameters: {info_to_update}")

        hoc_file = fileinput.input(
            files=os.path.join(self._neuron_workdir, "axon4pyfull.hoc"), inplace=1
        )
        for line in hoc_file:
            if any(line.startswith(matched_info := info) for info in info_to_update):
                _logger.debug(f"Matched: {matched_info}")
                _logger.debug(line)
                replacement_line = f"{matched_info}={parameters_dict[matched_info]}\n"
                line = replacement_line
                _logger.debug(line)
            print(line, end="")
        hoc_file.close()

    def modify_hoc_file(self, nRanvier, stepsPerMs, axon_morphology):
        axonDiam = axon_morphology["axon_diam"]
        if axonDiam >= 5.7:
            axoninter = (nRanvier - 1) * 6
        else:
            axoninter = (nRanvier - 1) * 3

        parameters_dict = {
            "axonnodes": nRanvier,
            "paranodes1": axon_morphology["n_para1"],
            "paranodes2": axon_morphology["n_para2"],
            "axoninter": axoninter,
            "axontotal": axon_morphology["n_segments"],
            "v_init": self._v_init,
            "fiberD": axonDiam,
            "paralength1": axon_morphology["para1_length"] * 1e3,
            "paralength2": axon_morphology["para2_length"] * 1e3,
            "nodelength": axon_morphology["ranvier_length"] * 1e3,
            "nodeD": axon_morphology["node_d"],
            "axonD": axon_morphology["axon_d"],
            "paraD1": axon_morphology["para1_d"],
            "paraD2": axon_morphology["para2_d"],
            "deltax": axon_morphology["node_step"] * 1e3,
            "nl": axon_morphology["lamellas"],
            "steps_per_ms": stepsPerMs,
        }
        self.paste_to_hoc(parameters_dict)

    def upsample_voltage(self, v_time_sol, axonDiam, axon_morphology):
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
        TODO refactor code
        """
        _logger.debug("Upsampling voltage")
        # let's interpolate voltage between node - center_l - center_r - node
        # assume 11 segments
        # n_segments_ds = ((n_segments_full - 1) / 11) * 3 +1

        n_segments_actual = axon_morphology["n_segments"]
        _logger.debug(
            f"Upsampling from {v_time_sol.shape[0]} to {n_segments_actual} segments"
        )

        v_time_sol_full = np.zeros((n_segments_actual, v_time_sol.shape[1]), float)

        if axonDiam >= 5.7:
            # fill out nodes first
            for k in np.arange(0, n_segments_actual, 11):
                z = int(k / 11) * 3
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # now two segments in between
            for k in np.arange(3, n_segments_actual, 11):
                z = int(k / 11) * 3 + 1
                v_time_sol_full[k, :] = v_time_sol[z, :]

            for k in np.arange(8, n_segments_actual, 11):
                z = int(k / 11) * 3 + 2
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # node -- -- intern -- -- -- -- intern -- -- node  ->  node-para1-para2-intern-intern-intern-intern-intern-intern-para2-para1-node
            internodal_length = (
                axon_morphology["node_step"]
                - axon_morphology["ranvier_length"]
                - (axon_morphology["para2_length"] + axon_morphology["para1_length"])
                * 2.0
            ) / 6.0
            dist_node_internode = (
                (axon_morphology["ranvier_length"] + internodal_length) * 0.5
                + axon_morphology["para1_length"]
                + axon_morphology["para2_length"]
            )
            ratio_1 = (
                (axon_morphology["ranvier_length"] + axon_morphology["para1_length"])
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (axon_morphology["para1_length"] + axon_morphology["para2_length"])
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
                for j in np.arange(0, n_segments_actual - 1, 11):
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
            for k in np.arange(0, n_segments_actual, 8):
                z = int(k / 8) * 2
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # now the center between nodes
            for k in np.arange(4, n_segments_actual, 8):
                z = int(k / 8) * 2 + 1
                v_time_sol_full[k, :] = v_time_sol[z, :]

            # node -- -- -- internodal -- -- -- node  ->  node-para1-para2-intern-intern-intern-para2-para1-node
            dist_node_internode = axon_morphology["node_step"] / 2.0
            internodal_length = (
                axon_morphology["node_step"]
                - axon_morphology["ranvier_length"]
                - (axon_morphology["para2_length"] + axon_morphology["para1_length"])
                * 2
            ) / 3.0
            ratio_1 = (
                (axon_morphology["ranvier_length"] + axon_morphology["para1_length"])
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (axon_morphology["para1_length"] + axon_morphology["para2_length"])
                * 0.5
                / dist_node_internode
            )
            ratio_3 = (
                ratio_2
                + (axon_morphology["para2_length"] + internodal_length)
                * 0.5
                / dist_node_internode
            )

            # now interpolate to the rest
            list_interp = [
                [1, 2, 3],
                [5, 6, 7],
            ]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, axon_morphology["n_segments"] - 1, 8):
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
    _hoc_file = "init_B5_extracellular.hoc"

    @property
    def resources_path(self):
        return "McNeal1976"

    def get_axon_morphology(
        self, axon_diam, axon_length=None, n_Ranvier=None, downsampled=False
    ) -> dict:
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

    def paste_to_hoc(self, parameters_dict):
        info_to_update = ["axonnodes", "v_init", "steps_per_ms"]

        if not set(info_to_update) == set(parameters_dict.keys()):
            raise ValueError("Need to provide all parameters: {info_to_update}")

        hoc_file = fileinput.input(
            files=os.path.join(self._neuron_workdir, "init_B5_extracellular.hoc"),
            inplace=1,
        )

        for line in hoc_file:
            if any(line.startswith(matched_info := info) for info in info_to_update):
                _logger.debug(f"Matched: {matched_info}")
                _logger.debug(line)
                replacement_line = f"{matched_info}={parameters_dict[matched_info]}\n"
                line = replacement_line
                _logger.debug(line)
            print(line, end="")
        hoc_file.close()

        NNODES_line = "NNODES ="
        axonnodes = parameters_dict["axonnodes"]
        NNODES_input = f"NNODES = {axonnodes}\n"

        hoc_file = fileinput.input(
            files=os.path.join(self._neuron_workdir, "axon5.hoc"), inplace=1
        )
        for line in hoc_file:
            if line.startswith(NNODES_line):
                line = NNODES_input
            print(line, end="")
        hoc_file.close()

        return True

    def modify_hoc_file(self, nRanvier, stepsPerMs, axon_morphology):
        parameters_dict = {
            "v_init": self._v_init,
            "axonnodes": nRanvier,
            "steps_per_ms": stepsPerMs,
        }
        self.paste_to_hoc(parameters_dict)
