"""Generate pam_results_summary.csv for the StimSets convergence study.

Reads the per-strategy overview CSVs (from evaluate_convergence.py)
and VCM_report.json files to produce a summary table with DOFs,
compute time, and per-pathway mean/max absolute error vs. best.
"""

import json
import os

import numpy as np
import pandas as pd

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
    "Results_PAM_default_meshsize",
    "Results_PAM_fine",
    "Results_PAM_material_refinement",
    "Results_PAM_fine_edge_refinement",
    "Results_PAM_edge_single_material_refinement",
    "Results_PAM_hp_refinement",
    "Results_PAM_hp_material_refinement",
]

best_directory = "Results_PAM_best"

# Load best (reference) results
best_df = pd.read_csv(f"{best_directory}_overview.csv")

# Get best DOF and timing
best_vcm = json.load(open(f"{best_directory}E1C1/VCM_report.json"))
best_dofs = best_vcm["DOF"]
best_time = sum(best_vcm["Timings"]["ComputeSolution"])
# Sum across all 8 contacts
for i in range(2, 9):
    vcm_path = f"{best_directory}E1C{i}/VCM_report.json"
    if os.path.isfile(vcm_path):
        vcm = json.load(open(vcm_path))
        best_time += sum(vcm["Timings"]["ComputeSolution"])

# Build summary
results_dict = {
    "roman": [],
    "study_name": [],
    "time": [],
    "dofs": [],
}
for pathway in pathways:
    results_dict[f"Pathway_status_{pathway}.json_activated"] = []
    results_dict[f"Pathway_status_{pathway}.json_activated_rel_error"] = []
    results_dict[f"Pathway_status_{pathway}.json_activated_max_error"] = []

for idx, result_dir in enumerate(result_directories):
    study_name = result_dir.replace("Results_PAM_", "")

    # DOF and timing from E1C1 VCM report
    vcm_path = f"{result_dir}E1C1/VCM_report.json"
    vcm = json.load(open(vcm_path))
    dofs = vcm["DOF"]
    time_total = sum(vcm["Timings"]["ComputeSolution"])
    for i in range(2, 9):
        vcm_path_i = f"{result_dir}E1C{i}/VCM_report.json"
        if os.path.isfile(vcm_path_i):
            vcm_i = json.load(open(vcm_path_i))
            time_total += sum(vcm_i["Timings"]["ComputeSolution"])

    results_dict["roman"].append(r"\rom{" f"{idx + 1}" "}")
    results_dict["study_name"].append(study_name)
    results_dict["dofs"].append(dofs)
    results_dict["time"].append(time_total)

    # Per-pathway activation and error (mean over 254 protocols)
    df = pd.read_csv(f"{result_dir}_overview.csv")
    for pathway in pathways:
        mean_activation = df[pathway].mean()
        abs_errors = np.abs(df[pathway].values - best_df[pathway].values)
        mean_error = abs_errors.mean()
        max_error = abs_errors.max()
        results_dict[f"Pathway_status_{pathway}.json_activated"].append(mean_activation)
        results_dict[f"Pathway_status_{pathway}.json_activated_rel_error"].append(
            mean_error
        )
        results_dict[f"Pathway_status_{pathway}.json_activated_max_error"].append(
            max_error
        )

summary = pd.DataFrame(results_dict)
summary["best_time"] = best_time
summary["best_dofs"] = best_dofs

for pathway in pathways:
    summary[f"Pathway_status_{pathway}.json_activated_best"] = best_df[pathway].mean()

summary.to_csv("pam_results_summary.csv", index=False)
print("Wrote pam_results_summary.csv")
print(summary[["roman", "study_name", "dofs", "time"]].to_string(index=False))
