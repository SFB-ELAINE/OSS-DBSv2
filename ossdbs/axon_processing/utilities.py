import json
import os

import numpy as np
from scipy.io import savemat


def create_leaddbs_outputs(
    output_path, Axon_Lead_DBS, connectome_name, scaling_index=None, pathway_name=None
):
    """Export axons with activation state in Lead-DBS supported format.

    Parameters
    ----------
    Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

    """
    mdic = {
        "fibers": Axon_Lead_DBS,
        "ea_fibformat": "1.0",
        "connectome_name": connectome_name,
    }  # For Lead-DBS .mat files

    if scaling_index is None:
        if pathway_name is None:
            savemat(os.path.join(output_path, "Axon_state.mat"), mdic)
        else:
            savemat(os.path.join(output_path, f"Axon_state_{pathway_name}.mat"), mdic)
    else:
        if pathway_name is None:
            savemat(os.path.join(output_path, f"Axon_state_{scaling_index}.mat"), mdic)
        else:
            savemat(
                os.path.join(
                    output_path, f"Axon_state_{pathway_name}_{scaling_index}.mat"
                ),
                mdic,
            )


def create_paraview_outputs(
    output_path, Axon_Lead_DBS, scaling_index=None, pathway_name=None
):
    """Export axons with activation state in Paraview supported format.

    Parameters
    ----------
    Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

    """
    if scaling_index is None:
        if pathway_name is None:
            np.savetxt(
                os.path.join(output_path, "Axon_state.csv"),
                Axon_Lead_DBS,
                delimiter=",",
                header="x-pt,y-pt,z-pt,idx,status",
            )
        else:
            np.savetxt(
                os.path.join(output_path, f"Axon_state_{pathway_name}.csv"),
                Axon_Lead_DBS,
                delimiter=",",
                header="x-pt,y-pt,z-pt,idx,status",
            )
    else:
        if pathway_name is None:
            np.savetxt(
                os.path.join(output_path, f"Axon_state_{scaling_index}.csv"),
                Axon_Lead_DBS,
                delimiter=",",
                header="x-pt,y-pt,z-pt,idx,status",
            )
        else:
            np.savetxt(
                os.path.join(
                    output_path, f"Axon_state_{pathway_name}_{scaling_index}.csv"
                ),
                Axon_Lead_DBS,
                delimiter=",",
                header="x-pt,y-pt,z-pt,idx,status",
            )


def store_axon_statuses(
    output_path,
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

    if scaling_index is None:
        if pathway_name is None:
            path_to_save = os.path.join(output_path, "Pathway_status.json")
        else:
            summary_dict["pathway_name"] = pathway_name
            path_to_save = os.path.join(
                output_path, f"Pathway_status_{pathway_name}.json"
            )
    else:
        summary_dict["scaling_index"] = str(scaling_index)
        if pathway_name is None:
            path_to_save = os.path.join(
                output_path, f"Pathway_status_{scaling_index}.json"
            )
        else:
            summary_dict["pathway_name"] = pathway_name
            path_to_save = os.path.join(
                output_path,
                f"Pathway_status_{pathway_name}_{scaling_index}.json",
            )
    with open(path_to_save, "w") as save_as_dict:
        json.dump(summary_dict, save_as_dict)
