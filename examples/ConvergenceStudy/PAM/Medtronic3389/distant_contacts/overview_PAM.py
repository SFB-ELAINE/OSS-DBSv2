import glob
import json
import os

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
info_str = f"{vcm_info['DOF']}"
for pathway in pathway_files:
    with open(os.path.join(best_PAM_result_directory, pathway)) as fp:
        pam_best[pathway] = json.load(fp)
    info_str += " "
    info_str += f"{pam_best[pathway]['percent_activated']}"

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


print("DOF", pathway_files)
print("Results_PAM_best")
print(info_str)
for result_dir in result_directories:
    print(result_dir)
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm_info = json.load(fp)
    info_str = f"{vcm_info['DOF']}"
    for pathway in pathway_files:
        with open(os.path.join(result_dir, pathway)) as fp:
            pam_info = json.load(fp)
        info_str += " "
        best_pct = pam_best[pathway]["percent_activated"]
        diff = abs(best_pct - pam_info["percent_activated"])
        info_str += f"{diff:.2f}"
    print(info_str)
