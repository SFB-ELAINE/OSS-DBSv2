import json
import os

import pandas as pd

stim_sets_combinations = 254
pathways = [
    "M1_cf_face_right",
    "M1_cf_lowerex_right",
    "M1_cf_upperex_right",
    "R_M1_hdp_face_right",
    "R_M1_hdp_lowerex_right",
    "R_M1_hdp_upperex_right",
    "cerebellothalamic_right",
    "gpe2stn_ass_right",
    "gpe2stn_sm_right",
    "medial_lemniscus_right",
]

result_directories = [
    "Results_PAM_best",
    "Results_PAM_default",
    "Results_PAM_fine",
    "Results_PAM_fine_edge_refinement",
    "Results_PAM_material_refinement",
    "Results_PAM_edge_single_material_refinement",
]

# get best result
for result_directory in result_directories:
    directory_results = {}
    for pathway in pathways:
        if pathway not in directory_results:
            directory_results[pathway] = []
        for stim_idx in range(stim_sets_combinations):
            with open(
                os.path.join(
                    result_directory, f"Pathway_status_{pathway}_{stim_idx}.json"
                )
            ) as fp:
                result_pathway = json.load(fp)
                directory_results[pathway].append(result_pathway["percent_activated"])
    df = pd.DataFrame(directory_results)
    df.to_csv(result_directory + "_overview.csv", index=False)
