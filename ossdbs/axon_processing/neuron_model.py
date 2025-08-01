import fileinput
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from importlib.resources import files
from typing import Optional

import h5py
import neuron
import numpy as np

from .axon_models import AxonMorphologyMcNeal1976, AxonMorphologyMRG2002
from .utilities import (
    create_leaddbs_outputs,
    create_paraview_outputs,
    store_axon_statuses,
)

_logger = logging.getLogger(__name__)


NEURON_DIR = "neuron_model"


class NeuronSimulator(ABC):
    """Interface to NEURON simulator."""

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
    # axon morphology class
    _ax_morph_class = None

    def __init__(
        self,
        pathways_dict: dict,
        output_path: str,
        scaling_vector: Optional[list] = None,
    ):
        self._neuron_executable = "nrnivmodl"
        # executable named different on Windows
        if sys.platform == "win32":
            self._neuron_executable = "mknrndll"
        self._pathways_dict = pathways_dict
        self._scaling_vector = None

        self._output_path = output_path
        # create output path if it does not exist
        if not os.path.isdir(output_path):
            _logger.info(f"Created directory {output_path}")
            os.mkdir(self.output_path)

        # directory, where all NEURON simulations will be conducted
        self._neuron_workdir = os.path.join(self._output_path, NEURON_DIR)
        # create local directory if it does not exist
        if not os.path.isdir(self._neuron_workdir):
            os.mkdir(self._neuron_workdir)

        # copy NEURON files and compile them
        _logger.info("Copy and compile NEURON files")
        self.copy_neuron_files()
        self.compile_neuron_files()

        # time-domain solution is loaded
        self._signal_dict = None

        # check if downsampled
        if "_DS" in pathways_dict["Axon_Model_Type"]:
            self._downsampled = True

    @property
    def hoc_file(self):
        """Name of hoc file."""
        if self._hoc_file is None:
            raise NotImplementedError("Need to implement default hoc file!")
        return self._hoc_file

    @property
    def axon_model(self):
        """Name of axon model."""
        return self._axon_model

    @property
    def signal_dict(self):
        """Information about stimulation signal."""
        return self._signal_dict

    @signal_dict.setter
    def signal_dict(self, value: dict):
        self._signal_dict = value

    @property
    def axon_model_type(self):
        """Type of axon model."""
        return self._pathways_dict["Axon_Model_Type"]

    @property
    def connectome_name(self):
        """Name of connectome."""
        return self._pathways_dict["connectome_name"]

    def get_axon_diam(self, idx: int):
        """Return the axon/fiber diameters."""
        return self._pathways_dict["axon_diams"][idx]

    def get_n_Ranvier(self, idx: int):
        """Get the number of nodes of Ranvier."""
        return self._pathways_dict["n_Ranvier"][idx]

    def get_N_seeded_neurons(self, idx: int):
        """Get the number of seeded neurons."""
        return self._pathways_dict["N_seeded_neurons"][idx]

    def get_N_orig_neurons(self, idx: int):
        """Get the number of original neurons."""
        return self._pathways_dict["N_orig_neurons"][idx]

    @property
    def downsampled(self) -> bool:
        """Return if model is downsampled."""
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

    def copy_neuron_files(self):
        """Copy files from template folder to local folder."""
        path_to_ossdbs = files("ossdbs")
        path_to_neuron_files = path_to_ossdbs.joinpath(
            "axon_processing", "neuron_templates", self.resources_path
        )
        # Directory can already exist
        shutil.copytree(path_to_neuron_files, self._neuron_workdir, dirs_exist_ok=True)

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

    def compile_neuron_files(self):
        """Compile a NEURON file."""
        _logger.info("Compile NEURON executable")
        subprocess.run(
            self.neuron_executable,
            # do not change! causes OSError
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            cwd=os.path.abspath(self._neuron_workdir),
        )
        _logger.info("Load mechanisms into environment")
        # TODO should be written in a safer way
        try:
            neuron.load_mechanisms(self._neuron_workdir)
        except RuntimeError:
            _logger.warning(
                "Mechanism was already loaded, continuing"
                " with loaded mechanism. "
                "Run a new Python instance to circumvent this error."
            )
            pass

    def load_solution(self, time_domain_h5_file: str):
        """Load solution from h5 file.

        Parameters
        ----------
        time_domain_h5_file: str
            Name of file holding time-domain solution

        Returns
        -------
        td_solution: h5 dataset
            Time-domain solution
        """
        td_solution = h5py.File(time_domain_h5_file, "r")

        return td_solution

    def load_unit_solutions(self, time_domain_h5_files: list):
        """Load solutions from h5 file for each contact-ground.

        Parameters
        ----------
        time_domain_h5_files: list
            names of files holding time-domain solution

        Returns
        -------
        td_unit_solutions: list
            unit time-domain solutions for each contact
        """
        td_unit_solutions = []  # list where we store datasets
        for solution_i in range(len(time_domain_h5_files)):
            td_solution = h5py.File(time_domain_h5_files[solution_i], "r")
            td_unit_solutions.append(td_solution)

        return td_unit_solutions

    def superimpose_unit_solutions(self, td_unit_solutions, scaling_vector: list):
        """Scale and superimpose solutions obtained for each contact-ground.
        This should be done across all datagroups and datasets.

        Parameters
        ----------
        td_unit_solutions: list
            unit time-domain solutions for each contact
        scaling_vector: list
            current scaling across contacts

        Returns
        -------
        td_solution: h5 dataset
            Superimposed and scaled time-domain solution
        """
        # very dumb way to get td_solution initialized
        td_solution = h5py.File(
            os.path.join(self.output_path, "combined_solution.h5"), mode="w"
        )
        for obj in td_unit_solutions[0].keys():
            td_unit_solutions[0].copy(obj, td_solution)

        pathways = list(td_unit_solutions[0].keys())
        pathways.remove("TimeSteps[s]")

        for pathway_idx, pathway_name in enumerate(pathways):
            N_neurons = self.get_N_seeded_neurons(pathway_idx)

            # Store Status once, assuming it's consistent
            status_dataset_first_solution = td_unit_solutions[0][pathway_name]["Status"]

            for neuron_index in range(N_neurons):
                if status_dataset_first_solution[neuron_index] == 0:
                    neuron_name = "axon" + str(neuron_index)

                    # Read all neuron potentials for this axon across all solutions into a list of arrays
                    neuron_potentials_all_solutions = []
                    for solution_i in range(len(td_unit_solutions)):
                        pathway_dataset = td_unit_solutions[solution_i][pathway_name]
                        neuron = pathway_dataset[neuron_name]
                        neuron_potentials_all_solutions.append(
                            np.array(neuron["Potential[V]"])
                        )

                    all_potentials_stacked = np.stack(
                        neuron_potentials_all_solutions, axis=0
                    )

                    scaling_vector_np = np.array(scaling_vector)

                    # Perform the vectorized scalar multiplication and sum
                    # (num_solutions, time_points) * (num_solutions,) -> (num_solutions, time_points)
                    # The multiplication broadcasts scaling_vector_np along axis=1
                    scaled_potentials = (
                        all_potentials_stacked
                        * scaling_vector_np[:, np.newaxis, np.newaxis]
                    )

                    # Sum along the solutions axis to get the superimposed result

                    td_solution[pathway_name]["axon" + str(neuron_index)][
                        "Potential[V]"
                    ][...] = np.sum(scaled_potentials, axis=0)

        return td_solution

    def process_pathways(
        self, td_solution, scaling: float = 1.0, scaling_index: Optional[int] = None
    ):
        """Go through all pathways and compute the activation.

        Parameters
        ----------
        td_solution: h5 dataset
            time-domain solution
        scaling: float
            scaling factor for the whole solution
        scaling_index: int
            index of the scaling factor or scaling vector

        """
        if td_solution is None:
            raise ValueError("Need to load time-domain solution from H5 file first.")
        pathways = list(td_solution.keys())
        pathways.remove("TimeSteps[s]")

        # signal parameters can be extracted from solution
        TimeSteps = np.array(td_solution["TimeSteps[s]"])

        self.signal_dict = {
            "time_step": np.round(
                1000.0 * (TimeSteps[1] - TimeSteps[0]), 6
            ),  # convert to ms
            "scaling": scaling,  # from Lead-DBS GUI
            "N_time_steps": TimeSteps.shape[0],
        }

        _logger.info("Going through pathways")
        for pathway_idx, pathway_name in enumerate(pathways):
            pathway_dataset = td_solution[pathway_name]
            self.check_pathway_activation(
                pathway_dataset, pathway_idx, pathway_name, scaling_index
            )

        # close H5 file
        td_solution.close()

    def get_v_ext(self, v_time_sol):
        """Convert potential computed by OSS-DBS to extracellular potential.

        Parameters
        ----------
        v_time_sol: numpy.ndarray
            NxM potential distribution (in V) in space (on the neuron)
            and time (DBS signal)

        Returns
        -------
        v_ext: numpy.ndarray
            NxM extracellular potential distribution (mV) in space (on the neuron)
            and time (DBS signal)

        Notes
        -----
        Scaling is taken into account!
        """
        scaling = self.signal_dict["scaling"]
        v_ext = np.zeros_like(v_time_sol, float)
        if self._scaling_vector is None:
            return v_time_sol * 1000.0 * scaling  # convert to mV

        # case with contact-wise scaling
        for scaling_comp in self._scaling_vector:
            v_ext = v_ext + v_time_sol * 1000.0 * scaling * scaling_comp
        return v_ext

    def get_axon_status_multiprocessing(self, neuron_index, v_time_sol, output):
        """Probe action potential at a neuron (axon).

        Parameters
        ----------
        neuron_index: int
            index of neuron in the pathway starting from 0
        v_time_sol: np.ndarray
            potential distribution (in V) in space (on the neuron)
            and time (DBS signal)
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
        _logger.debug(f"Load file: {os.path.join(self._neuron_workdir, self.hoc_file)}")
        neuron.h(
            f'{{load_file("{os.path.join(self._neuron_workdir, self.hoc_file)}")}}'
        )
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

    # ruff: noqa: C901
    def check_pathway_activation(
        self, pathway_dataset, pathway_idx, pathway_name=None, scaling_index=None
    ):
        """Parallelized probing of action potentials at all stimulated neurons
           (axons) described by supplied pathway datagroups.

        Parameters
        ----------
        pathway_dataset: h5 group
            contains datasets that describe geometries for all neuron (axons)
        pre_status: numpy.ndarray
            vector of initial neuron (axon) statuses
            (0: available for probing, -1: intersected with electrode, -2: inside CSF)
        pathway_name: str, optional
            Name of pathway
        pathway_idx: int
            Index of pathway
        scaling_index: int, optional
            Index of scaling (needed for superposition)

        Note
        ----------
        Creates 'Axon_states_*' .mat and .csv files for visualization
        in Lead-DBS and Paraview, respectively.
        Also stores summary statistics in 'Pathway_status_*.json'
        """
        # use half of CPUs
        N_proc = mp.cpu_count() / 2

        # get parameters
        N_neurons = self.get_N_seeded_neurons(pathway_idx)
        axon_diam = self.get_axon_diam(pathway_idx)
        n_Ranvier = self.get_n_Ranvier(pathway_idx)
        orig_N_neurons = self.get_N_orig_neurons(pathway_idx)

        axon_morphology = self._ax_morph_class()
        axon_morphology.update_axon_morphology(axon_diam, n_Ranvier=n_Ranvier)
        if self.downsampled:
            # check actual number of n_segments in case downsampled
            n_segments_actual = axon_morphology.get_n_segments(
                downsampled=self.downsampled
            )
            _logger.debug(f"n_segments_actual: {n_segments_actual}")
        else:
            n_segments_actual = axon_morphology.n_segments
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
                    neuron_index * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    :3,
                ] = np.array(neuron["Points[mm]"])

                # add index
                Axon_Lead_DBS[
                    neuron_index * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    3,
                ] = neuron_index + 1  # because Lead-DBS numbering starts from 1

                # check which neurons were flagged with CSF and electrode intersection
                # skip probing of those
                if pre_status[neuron_index] != 0:
                    Axon_Lead_DBS[
                        neuron_index * n_segments_actual : (neuron_index + 1)
                        * n_segments_actual,
                        4,
                    ] = pre_status[neuron_index]
                    neuron_index += 1
                    continue

                neuron_time_sol = np.array(neuron["Potential[V]"])
                # upsample, TODO maybe move inside mp.Process
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
                    neuron_index * n_segments_actual : (neuron_index + 1)
                    * n_segments_actual,
                    4,
                ] = 1
            elif neuron_index in List_of_not_activated:
                Axon_Lead_DBS[
                    neuron_index * n_segments_actual : (neuron_index + 1)
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

        if scaling_index is None:
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
    """MRG2002 NEURON model."""

    _axon_model = "MRG2002"
    # needs extra compilation steps
    _extra_initialization = True
    # different initial resting potential
    _v_init = -80.0
    _hoc_file = "axon4pyfull.hoc"

    _ax_morph_class = AxonMorphologyMRG2002

    @property
    def resources_path(self):
        """Link to template NEURON files."""
        return "MRG2002"

    def paste_to_hoc(self, parameters_dict: dict) -> None:
        """Update the hoc file with parameters."""
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
        """Update parameters in the hoc file."""
        fiber_diam = axon_morphology.fiber_diam
        if fiber_diam >= 5.7:
            axoninter = (nRanvier - 1) * 6
        else:
            axoninter = (nRanvier - 1) * 3

        parameters_dict = {
            "axonnodes": nRanvier,
            "paranodes1": axon_morphology.n_para1,
            "paranodes2": axon_morphology.n_para2,
            "axoninter": axoninter,
            "axontotal": axon_morphology.n_segments,
            "v_init": self._v_init,
            "fiberD": fiber_diam,
            "paralength1": axon_morphology.para1_length * 1e3,
            "paralength2": axon_morphology.para2_length * 1e3,
            "nodelength": axon_morphology.ranvier_length * 1e3,
            "nodeD": axon_morphology.node_d,
            "axonD": axon_morphology.axon_d,
            "paraD1": axon_morphology.para1_d,
            "paraD2": axon_morphology.para2_d,
            "deltax": axon_morphology.node_step * 1e3,
            "nl": axon_morphology.lamellas,
            "steps_per_ms": stepsPerMs,
        }
        self.paste_to_hoc(parameters_dict)

    def upsample_voltage(self, v_time_sol, axonDiam, axon_morphology):
        """Upsample potential distribution for downsampled neurons by interpolating.

        Parameters
        ----------
        v_time_sol: numpy.ndarray
            LxM potential distribution (in V) in space (on the neuron)
            and time (DBS signal)
        axonDiam: float
            Diameter of axon
        axon_morphology: AxonMorphology
            Model of axon morphology

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

        n_segments_actual = axon_morphology.n_segments
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

            # node -- -- intern -- -- -- -- intern -- -- node ->
            # ruff: noqa: E501
            # node-para1-para2-intern-intern-intern-intern-intern-intern-para2-para1-node
            internodal_length = (
                axon_morphology.node_step
                - axon_morphology.ranvier_length
                - (axon_morphology.para2_length + axon_morphology.para1_length) * 2.0
            ) / 6.0
            dist_node_internode = (
                (axon_morphology.ranvier_length + internodal_length) * 0.5
                + axon_morphology.para1_length
                + axon_morphology.para2_length
            )
            ratio_1 = (
                (axon_morphology.ranvier_length + axon_morphology.para1_length)
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (axon_morphology.para1_length + axon_morphology.para2_length)
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
            dist_node_internode = axon_morphology.node_step / 2.0
            internodal_length = (
                axon_morphology.node_step
                - axon_morphology.ranvier_length
                - (axon_morphology.para2_length + axon_morphology.para1_length) * 2
            ) / 3.0
            ratio_1 = (
                (axon_morphology.ranvier_length + axon_morphology.para1_length)
                * 0.5
                / dist_node_internode
            )
            ratio_2 = (
                ratio_1
                + (axon_morphology.para1_length + axon_morphology.para2_length)
                * 0.5
                / dist_node_internode
            )
            ratio_3 = (
                ratio_2
                + (axon_morphology.para2_length + internodal_length)
                * 0.5
                / dist_node_internode
            )

            # now interpolate to the rest
            list_interp = [
                [1, 2, 3],
                [5, 6, 7],
            ]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, axon_morphology.n_segments - 1, 8):
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
    """McNeal1976 NEURON model."""

    _axon_model = "McNeal1976"
    _hoc_file = "init_B5_extracellular.hoc"
    _ax_morph_class = AxonMorphologyMcNeal1976

    @property
    def resources_path(self):
        """Link to template NEURON files."""
        return "McNeal1976"

    def paste_to_hoc(self, parameters_dict):
        """Update the hoc file with parameters."""
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
        """Update parameters in the hoc file."""
        parameters_dict = {
            "v_init": self._v_init,
            "axonnodes": nRanvier,
            "steps_per_ms": stepsPerMs,
        }
        self.paste_to_hoc(parameters_dict)
