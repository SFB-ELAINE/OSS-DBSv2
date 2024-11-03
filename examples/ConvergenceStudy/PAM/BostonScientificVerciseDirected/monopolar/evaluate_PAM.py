import glob
import os
from pprint import pprint

import pandas as pd

from ossdbs.axon_processing import compare_pathways

best_PAM_result_directory = "Results_PAM_best"
pam_best_files = glob.glob(os.path.join(best_PAM_result_directory, "Axon_state*.csv"))

result_directories = [
    "Results_PAM_default",
    "Results_PAM_default_meshsize",
    "Results_PAM_fine",
    "Results_PAM_very_fine",
    "Results_PAM_edge_refinement",
    "Results_PAM_fine_edge_refinement",
    "Results_PAM_very_fine_edge_refinement",
    "Results_PAM_edge_voxel_refinement",
    "Results_PAM_edge_meshsize",
    "Results_PAM_material_refinement",
    "Results_PAM_edge_single_material_refinement",
    "Results_PAM_edge_double_material_refinement",
]


def compare_pathway_files(pam_files: list, pam_best_directory: str):
    """Compare pathway files from list. Load corresponding file
    from directory where benchmark result is stored.
    """
    comparison = {}
    for pam_file in pam_files:
        axon_file_name = pam_file.split(os.sep)[-1]
        pam_best_file = os.path.join(pam_best_directory, axon_file_name)
        pam_data = pd.read_csv(pam_file)
        pam_best_data = pd.read_csv(pam_best_file)
        comparison[axon_file_name] = compare_pathways(pam_data, pam_best_data)
        print(axon_file_name)
        pprint(comparison[axon_file_name])
    return comparison


for result_dir in result_directories:
    pam_files = glob.glob(os.path.join(result_dir, "Axon_state*.csv"))

    print("###############\n\n")
    print(
        result_dir.replace("Results_PAM_", ""),
    )
    comparison = compare_pathway_files(pam_files, best_PAM_result_directory)
    print("###############\n\n")
