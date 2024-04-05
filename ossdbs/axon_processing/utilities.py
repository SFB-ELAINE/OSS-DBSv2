# Copyright 2024 Konstantin Butenko, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import os

import numpy as np
from dipy.tracking.streamline import set_number_of_points
from nibabel.streamlines import ArraySequence
from scipy.io import savemat


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def convert_fibers_to_streamlines(fibers):
    """Convert Lead-DBS fibers to Nibabel streamlines.

    Parameters
    ----------
    fibers,: fiber descriptions (Lead-DBS format)

    Returns
    -------
    list, streamlines stored as ArraySequence(), i.e. list that describes each fiber in a sublist

    """
    streamlines = ArraySequence()

    # yes, indexing starts with one in those .mat files
    N_streamlines = int(fibers[3, :].max())

    k = 0
    i_previous = 0
    for i in range(N_streamlines):
        loc_counter = 0
        while (i + 1) == fibers[
            3, k
        ]:  # this is not optimal, you need to extract a pack by np.count?
            k += 1
            loc_counter += 1
            if k == fibers[3, :].shape[0]:
                break

        stream_line = fibers[:3, i_previous : i_previous + loc_counter].T
        i_previous = k
        streamlines.append(stream_line)

    return streamlines


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


def resample_streamline_for_Ranvier(streamline_array, estim_axon_length, n_Ranvier):
    """Resample arbitrary streamline to equidistantly spaced nodes of Ranvier.

    Parameters
    ----------
    streamline_array:: list
        Each item in the list is an array with 3D coordinates of a streamline.
    estim_axon_length: float
        estimated length of axon for this streamline
    n_Ranvier: int
        Number of nodes of Ranviers at the seeded axon

    Returns
    -------
    list, streamline resampled to nodes of Ranvier
    """
    # after this index we cut the streamline (do not mix up with truncation to the actual axon!)
    cut_index, cummulated_length = index_for_length(streamline_array, estim_axon_length)

    # Don't mix up sums and positions. +1 for the last Ranvier node, +1 for the sum, +1 for index
    streamline_array_Ranvier = np.zeros((cut_index + 1 + 1 + 1, 3), float)
    last_segment_length = estim_axon_length - cummulated_length

    # adjust the last segment to match the estimated axon length exactly
    v = streamline_array[cut_index + 1 + 1, :] - streamline_array[cut_index + 1, :]
    v_hat = v / (v**2).sum() ** 0.5
    streamline_array_Ranvier[: cut_index + 1 + 1, :] = streamline_array[
        : cut_index + 1 + 1, :
    ]
    streamline_array_Ranvier[cut_index + 1 + 1, :] = (
        last_segment_length * v_hat + streamline_array[cut_index + 1, :]
    )

    streamline_resampled = set_number_of_points(
        streamline_array_Ranvier, nb_points=n_Ranvier
    )

    return streamline_resampled


def index_for_length(xyz, req_length, along=True):
    """Find streamline truncation point to match the given length.

    Parameters
    ----------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    req_length: float
        required length after truncation
    along : bool, optional
       If True, return array giving cumulative length along track,
       otherwise (default) return scalar giving total length.

    Returns
    -------
    int
        index of the streamline truncation point
    float
        length of the truncated streamline
    """
    xyz = np.asarray(xyz)
    if xyz.shape[0] < 2:
        if along:
            return np.array([0])
        return 0

    dists = np.sqrt((np.diff(xyz, axis=0) ** 2).sum(axis=1))

    if along:
        cummulated_lengths = np.cumsum(dists)
        idx, value = find_nearest(cummulated_lengths, req_length)
        if value > req_length:
            idx = idx - 1

        return idx, cummulated_lengths[idx]

    return idx, cummulated_lengths[idx]
