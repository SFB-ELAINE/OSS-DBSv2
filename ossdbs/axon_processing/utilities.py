# Copyright 2024 Konstantin Butenko, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import logging
import os
from typing import Optional

import numpy as np
import pandas as pd
from dipy.tracking.metrics import length as dipy_length
from dipy.tracking.streamline import set_number_of_points
from nibabel.streamlines import ArraySequence
from scipy import spatial
from scipy.io import savemat

_logger = logging.getLogger(__name__)


def compare_pathways(pam_df: pd.DataFrame, pam_df_to_compare: pd.DataFrame) -> dict:
    """Use data loaded from CSV files to compare axon activation."""
    try:
        import sklearn.metrics as metrics
    except ImportError as exc:
        raise ImportError("Please install scikit-learn to use their metrics.") from exc

    status = pam_df["status"].to_numpy()
    status_to_compare = pam_df_to_compare["status"].to_numpy()
    if np.all(np.isclose(status, 0)) and np.all(np.isclose(status_to_compare, 0)):
        _logger.warning("No axon activation found, returning None")
        return None
    if not status.shape == status_to_compare.shape:
        raise ValueError(
            "The provided DataFrames do not contain the same number of axons."
        )
    confusion_matrix = metrics.confusion_matrix(
        status, status_to_compare, labels=[-1, 0, 1]
    )
    precision, recall, f1, support = metrics.precision_recall_fscore_support(
        status, status_to_compare
    )
    return {
        "confusion matrix": confusion_matrix,
        "recall:": recall,
        "precision": precision,
        "f1": f1,
        "support": support,
    }


def find_nearest(array: np.ndarray, value: float):
    """Find index and value closest to value."""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def convert_fibers_to_streamlines(fibers: np.ndarray):
    """Convert Lead-DBS fibers to Nibabel streamlines.

    Parameters
    ----------
    fibers: np.ndarray
        fiber descriptions (Lead-DBS format)

    Returns
    -------
    list
        streamlines stored as ArraySequence(),
        i.e. list that describes each fiber in a sublist

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
    output_path: str,
    Axon_Lead_DBS: np.ndarray,
    connectome_name: str,
    scaling_index: Optional[int] = None,
    pathway_name: Optional[str] = None,
):
    """Export axons with activation state in Lead-DBS supported format.

    Parameters
    ----------
    output_path: str
        Path to save file
    Axon_Lead_DBS: numpy.ndarray
        NxM array with geometry, index and activation status of neurons
        (equivalent of connectome.fibers format in Lead-DBS)
    connectome_name: str
        Name of connectome_name
    scaling_index: int
        Index of scaling (used for superposition)
    pathway_name: str
        Name of pathway
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
    output_path: str,
    Axon_Lead_DBS: np.ndarray,
    scaling_index: Optional[int] = None,
    pathway_name: Optional[str] = None,
):
    """Export axons with activation state in Paraview supported format.

    Parameters
    ----------
    output_path: str
        Path to save file
    Axon_Lead_DBS: numpy.ndarray
        NxM array with geometry, index and activation status of neurons
        (equivalent of connectome.fibers format in Lead-DBS)
    connectome_name: str
        Name of connectome_name
    scaling_index: int
        Index of scaling (used for superposition)
    pathway_name: str
        Name of pathway

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
    output_path: str,
    percent_activated: float,
    percent_damaged: float,
    percent_csf: float,
    scaling_index: Optional[int] = None,
    pathway_name: Optional[bool] = None,
):
    """Store PAM results.

    Parameters
    ----------
    percent_activated: float
        percentage of the original(!) number of neurons that are activated
    percent_damaged: float
        percent of the original(!) number of neurons that are 'damaged'
    percent_csf: float
        percentage of the original(!) number of neurons that intersect with CSF
    output_path: str
        Path to save file
    scaling_index: int
        Index of scaling (used for superposition)
    pathway_name: str
        Name of pathway


    Notse
    -----
    For activation state of particular neuron see 'fiberActivation*' files
    as those restore original(!) indices as in the connectome.
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
    # after this index we cut the streamline
    # (do not mix up with truncation to the actual axon!)
    cut_index, cummulated_length = index_for_length(streamline_array, estim_axon_length)

    # Don't mix up sums and positions.
    # +1 for the last Ranvier node, +1 for the sum, +1 for index
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


def resample_fibers_to_Ranviers(streamlines: list, node_step: int, n_Ranvier: int):
    """Get streamlines resampled by nodes of Ranvier for a specific axonal morphology.

    Parameters
    ----------
    streamlines: list
        arbitrary sampled streamlines, stored as ArraySequence()
    node_step: int
        Length from a node of Ranvier to the next
    n_Ranvier: int
        Number of nodes of Ranvier

    Returns
    -------
    list, resampled streamlines, stored as ArraySequence()

    """
    # resampling to nodes of Ranvier for arbitrary fiber length
    lengths_streamlines_filtered = list(map(dipy_length, streamlines))
    streamlines_resampled = ArraySequence()

    excluded_streamlines = []
    # total_points = 0
    for streamline_index in range(len(lengths_streamlines_filtered)):
        n_Ranvier_this_axon = int(
            lengths_streamlines_filtered[streamline_index] / node_step
        )
        streamline_resampled = resample_streamline_for_Ranvier(
            streamlines[streamline_index],
            n_Ranvier_this_axon * node_step,
            n_Ranvier_this_axon,
        )
        if len(streamline_resampled) < n_Ranvier:
            _logger.info(f"Streamline {streamline_index} is too short")
            excluded_streamlines.append(streamline_index)
        else:
            streamlines_resampled.append(streamline_resampled)
            # total_points = total_points + len(streamline_resampled)

    return streamlines_resampled, excluded_streamlines


def normalized(vector: np.ndarray, axis: int = -1, order: int = 2):
    """Get L2 norm of a vector."""
    l2 = np.atleast_1d(np.linalg.norm(vector, order, axis))
    l2[l2 == 0] = 1
    return vector / np.expand_dims(l2, axis)


# ruff: noqa: C901
def place_axons_on_streamlines(
    streamlines_resampled: list, n_Ranvier: int, centering_coordinates: list
):
    """Allocate axons on streamlines at seeding points given by centering_coordinates.

    Parameters
    ----------
    streamlines_resampled: list
        streamlines sampled by nodes of Ranvier, stored as ArraySequence()
    n_Ranvier: int
        Number of nodes of Ranvier
    centering_coordinates: list of lists
        3-D coordinates used to center axons on fibers (e.g. active contacts)

    Returns
    -------
    list, axons (truncated streamlines), stored as ArraySequence()

    """
    axons_ROI_centered = ArraySequence()

    for inx_axn in range(len(streamlines_resampled)):
        single_streamline_ROI_centered = np.zeros((n_Ranvier, 3), float)

        A = streamlines_resampled[inx_axn]
        distance_list = []
        index_list = []
        for j in range(len(centering_coordinates)):
            distance, index = spatial.KDTree(A).query(
                centering_coordinates[j]
            )  # distance is a local index of closest node of Ranvier on the axon
            distance_list.append(distance)
            index_list.append(index)

        index = index_list[
            distance_list.index(min(distance_list))
        ]  # index of the closest point as assigned as index

        loc_index = 0
        # choose where to start seeding the axon
        if index < int(n_Ranvier / 2):
            # axon---fiber---fiber---fiber---fiber---#
            for i in range(0, int(n_Ranvier)):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index += 1
        elif index + int(n_Ranvier / 2) + 1 > A.shape[0]:
            # fiber---fiber---fiber---fiber---axon---#
            for i in range(A.shape[0] - n_Ranvier, A.shape[0]):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index += 1
        else:
            # ---fiber---fiber---axon---fiber---fiber---#
            if n_Ranvier % 2 == 0:
                for i in range(
                    index - int(n_Ranvier / 2),
                    index + int(n_Ranvier / 2),
                ):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index += 1
            else:
                for i in range(
                    index - int(n_Ranvier / 2),
                    index + int(n_Ranvier / 2) + 1,
                ):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index += 1

        axons_ROI_centered.append(single_streamline_ROI_centered)

    if len(axons_ROI_centered) != len(streamlines_resampled):
        raise RuntimeError("Failed to sample some axons!")

    return axons_ROI_centered
