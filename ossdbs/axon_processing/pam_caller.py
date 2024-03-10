import json
import logging
import os
import subprocess
import sys

import h5py
import numpy as np

from .axon_files.neuron_simulation import NeuronStimulation

_logger = logging.getLogger(__name__)


def call_pam(
    neuron_folder: str,
    folder_to_save: str,
    points_h5_file: str,
    pathways_params_file: str,
    scaling: float = 1.0,
    scaling_index=None,
):
    """Call to probe action potentials for a given time domain solution.

    Parameters
    ----------
    neuron_folder: str
        path to folder where NEURON models stored
    folder_to_save: str
        path to folder where results are stored. Lead-DBS expects <stim_folder>/Results_<hemis>
    points_h5_file: str
        path to .h5 containing the time domain solution for the pathways (point model)
    pathways_params_file: str
        path to .json containing parameters for the pathways
    scaling: float
        optional, scaling factor for the whole solution (different from scaling_vector)
    scaling_index: int
        optional, index of the scaling factor or scaling vector
    """
    # load files
    with open(pathways_params_file) as fp:
        pathways_dict = json.load(fp)

    neuron_executable = "nrnivmodl"
    # executable named different on Windows
    if sys.platform == "win32":
        neuron_executable = "mknrndll"

    # get to the right NEURON folder and compile
    if pathways_dict["Axon_Model_Type"] == "McNeal1976":
        _logger.info("Calling McNeal1976 model.")
        os.chdir(neuron_folder + "/McNeal1976")
        with open(os.devnull, "w") as FNULL:
            subprocess.call(
                neuron_executable, shell=True, stdout=FNULL, stderr=subprocess.STDOUT
            )
    elif pathways_dict["Axon_Model_Type"] in ["MRG2002", "MRG2002_DS"]:
        _logger.info("Calling MRG2002 model.")
        os.chdir(neuron_folder)
        with open(os.devnull, "w") as FNULL:
            subprocess.call(
                "nocmodl axnode.mod", shell=True, stdout=FNULL, stderr=subprocess.STDOUT
            )  # might not work with remote hard drives
            # TODO test if works
            subprocess.call(
                neuron_executable, shell=True, stdout=FNULL, stderr=subprocess.STDOUT
            )

    _logger.info("NEURON models compiled.")

    # load solution
    hf = h5py.File(points_h5_file, "r")
    pathways = list(hf.keys())
    pathways.remove("TimeSteps[s]")

    # signal parameters can be extracted from solution
    TimeSteps = np.array(hf["TimeSteps[s]"])
    signal_dict = {
        "time_step": np.round(1000.0 * (TimeSteps[1] - TimeSteps[0]), 6),  # in ms
        "scaling": scaling,  # from GUI
        "N_time_steps": TimeSteps.shape[0],
    }

    pathway_idx = 0
    _logger.info("Going through pathways")
    for pathway_name in pathways:
        pathway_dataset = hf[pathway_name]
        pathway_dict = {
            "pathway_name": pathway_name,
            "Axon_Model_Type": pathways_dict["Axon_Model_Type"],
            "axon_diam": pathways_dict["axon_diams"][pathway_idx],
            "n_Ranvier": pathways_dict["n_Ranvier"][pathway_idx],
            "N_seeded_neurons": pathways_dict["N_seeded_neurons"][pathway_idx],
            "N_orig_neurons": pathways_dict["N_orig_neurons"][pathway_idx],
            "connectome_name": pathways_dict["connectome_name"],
        }
        pathwayNEURON = NeuronStimulation(
            pathway_dict, signal_dict, folder_to_save, None, scaling_index
        )
        pathwayNEURON.check_pathway_activation(pathway_dataset)
        pathway_idx += 1
    hf.close()
