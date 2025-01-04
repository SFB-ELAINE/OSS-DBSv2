import glob
import json
import os

import pandas as pd

best_PAM_result_directory = "Results_PAM_best"
pathway_files = glob.glob(
    os.path.join(best_PAM_result_directory, "Pathway_status*.json")
)
pathway_files = [
    p.replace(best_PAM_result_directory + os.sep, "") for p in pathway_files
]
pam_best = {}
with open(os.path.join(best_PAM_result_directory, "VCM_report.json")) as fp:
    vcm_info = json.load(fp)
best_dofs = vcm_info["DOF"]
best_time = 0.0
for t in vcm_info["Timings"]["ComputeSolution"]:
    best_time += t
info_str = f'{vcm_info["DOF"]}'
results_dict = {}
for pathway in pathway_files:
    with open(os.path.join(best_PAM_result_directory, pathway)) as fp:
        pam_best[pathway] = json.load(fp)
    info_str += " "
    info_str += f'{pam_best[pathway]["percent_activated"]}'
    results_dict[f"{pathway}_activated"] = []
    results_dict[f"{pathway}_activated_rel_error"] = []

categories = ["roman", "study_name", "time", "dofs"]
for category in categories:
    results_dict[category] = []


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
    "Results_PAM_fine_edge_single_material_refinement",
    "Results_PAM_fine_edge_double_material_refinement",
]


print("DOF", pathway_files)
print("Results_PAM_best")
print(info_str)
for idx, result_dir in enumerate(result_directories):
    info_str = result_dir.replace("Results_PAM_", "")
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm_info = json.load(fp)
    info_str += " & "
    info_str += f'{vcm_info["DOF"]}'
    for pathway in pathway_files:
        with open(os.path.join(result_dir, pathway)) as fp:
            pam_info = json.load(fp)
        difference_in_activation = abs(
            pam_best[pathway]["percent_activated"] - pam_info["percent_activated"]
        )
        info_str += " & "
        info_str += f'{pam_info["percent_activated"]:.2f}'
        info_str += " & "
        info_str += f"{difference_in_activation:.2f}"
        results_dict[f"{pathway}_activated"].append(pam_info["percent_activated"])
        results_dict[f"{pathway}_activated_rel_error"].append(
            abs(pam_best[pathway]["percent_activated"] - pam_info["percent_activated"])
        )
    info_str += "\\\\"
    print(info_str)
    dofs = vcm_info["DOF"]
    time = 0.0
    for t in vcm_info["Timings"]["ComputeSolution"]:
        time += t
    study_name = result_dir.replace("Results_PAM_", "")
    results_dict["roman"].append(r"\rom{" f"{idx + 1}" "}")
    results_dict["dofs"].append(dofs)
    results_dict["time"].append(time)
    results_dict["study_name"].append(study_name)
df = pd.DataFrame(results_dict)
df["best_time"] = best_time
df["best_dofs"] = best_dofs

for pathway in pathway_files:
    df[f"{pathway}_activated_best"] = pam_best[pathway]["percent_activated"]
df.to_csv("pam_results_summary.csv", index=False)
