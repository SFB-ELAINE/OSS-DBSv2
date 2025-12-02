import numpy as np
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
    "Results_PAM_default",
    "Results_PAM_fine",
    "Results_PAM_fine_edge_refinement",
    "Results_PAM_material_refinement",
    "Results_PAM_edge_single_material_refinement",
]

best_result = pd.read_csv("Results_PAM_best_overview.csv")
# get best result
mean_activations = []
for result_directory in result_directories:
    print("++++++++++++++")
    print(result_directory)
    print("--------------")
    print("Pathway, Mean activation, Mean difference, Max. difference")
    df = pd.read_csv(result_directory + "_overview.csv")
    for pathway in pathways:
        best_activation = best_result[pathway].to_numpy()
        current_activation = df[pathway].to_numpy()
        difference = np.abs(best_activation - current_activation)
        mean_activation = np.mean(current_activation)
        std_activation = np.std(current_activation)
        mean_diff = np.mean(difference)
        std_diff = np.std(difference)
        max_diff = np.max(difference)
        print(f"{pathway}, {mean_activation:.2f}, {mean_diff:.2f}, {max_diff:.2f}")
    print("++++++++++++++")
