import json
import os

import numpy as np
import pandas as pd

from ossdbs import VTAImage

vta_voxel_best = VTAImage("Results_VTA_best/VTA_solution_Lattice.nii")
best_volume = vta_voxel_best.get_vta_volume()
with open(os.path.join("Results_VTA_best", "VCM_report.json")) as fp:
    vcm_info = json.load(fp)

best_dofs = vcm_info["DOF"]
best_time = vcm_info["Timings"]["ComputeSolution"][0]
with open(os.path.join("Results_VTA_best", "Lattice.json")) as fp:
    lattice_info = json.load(fp)
best_ngs_vta_volume = float(lattice_info["volume"])
impedance_data = pd.read_csv(os.path.join("Results_VTA_best", "impedance.csv"))
best_imp = impedance_data["real"].iloc[0]

print(
    "Best values: ",
    np.round(best_ngs_vta_volume, 2),
    np.round(best_volume, 2),
    np.round(best_imp, 2),
)

result_directories = [
    "Results_VTA_default",
    "Results_VTA_fine",
    "Results_VTA_very_fine",
    "Results_VTA_edge_refinement",
    "Results_VTA_fine_edge_refinement",
    "Results_VTA_very_fine_edge_refinement",
    "Results_VTA_edge_voxel_refinement",
    "Results_VTA_material_refinement",
    "Results_VTA_edge_single_material_refinement",
    "Results_VTA_edge_double_material_refinement",
    "Results_VTA_fine_edge_single_material_refinement",
    "Results_VTA_fine_edge_double_material_refinement",
]

print(
    "Directory &  DOF & time & Impedance & Rel. error & NGS VTA "
    "& Rel. error & VTA & Rel.error & Dice"
)
results_dict = {}
categories = [
    "roman",
    "study_name",
    "time",
    "dofs",
    "dice",
    "vta_voxel_volume",
    "ngs_vta_volume",
    "imp",
    "imp_rel_error",
    "ngs_vta_volume_rel_error",
    "vta_voxel_volume_rel_err",
]
for category in categories:
    results_dict[category] = []
for idx, result_dir in enumerate(result_directories):
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm_info = json.load(fp)
    with open(os.path.join(result_dir, "Lattice.json")) as fp:
        lattice_info = json.load(fp)
    impedance_data = pd.read_csv(os.path.join(result_dir, "impedance.csv"))
    imp = impedance_data["real"].iloc[0]
    imp_rel_error = 100.0 * abs(imp - best_imp) / best_imp
    ngs_vta_volume = float(lattice_info["volume"])
    vta_voxel = VTAImage(os.path.join(result_dir, "VTA_solution_Lattice.nii"))
    dice = vta_voxel_best.compute_dice_coefficent(vta_voxel)
    volume = vta_voxel.get_vta_volume()
    dofs = vcm_info["DOF"]
    time = vcm_info["Timings"]["ComputeSolution"][0]
    study_name = result_dir.replace("Results_VTA_", "")
    ngs_vta_volume_rel_error = (
        100.0 * abs(ngs_vta_volume - best_ngs_vta_volume) / best_ngs_vta_volume
    )
    volume_rel_error = 100.0 * abs(volume - best_volume) / best_volume
    print(
        study_name,
        " & ",
        dofs,
        " & ",
        np.round(time, 2),
        " & ",
        np.round(imp, 2),
        " & ",
        np.round(imp_rel_error, 2),
        " & ",
        np.round(ngs_vta_volume, 2),
        " & ",
        np.round(ngs_vta_volume_rel_error, 2),
        " & ",
        np.round(volume, 3),
        " & ",
        np.round(volume_rel_error, 2),
        " & ",
        np.round(dice, 4),
        " \\\\",
    )
    results_dict["roman"].append(r"\rom{" f"{idx + 1}" "}")
    results_dict["dofs"].append(dofs)
    results_dict["time"].append(time)
    results_dict["study_name"].append(study_name)
    results_dict["imp"].append(imp)
    results_dict["imp_rel_error"].append(imp_rel_error)
    results_dict["ngs_vta_volume"].append(ngs_vta_volume)
    results_dict["ngs_vta_volume_rel_error"].append(ngs_vta_volume_rel_error)
    results_dict["vta_voxel_volume"].append(volume)
    results_dict["vta_voxel_volume_rel_err"].append(volume_rel_error)
    results_dict["dice"].append(dice)

df = pd.DataFrame(results_dict)
df["best_impedance"] = best_imp
df["best_time"] = best_time
df["best_dofs"] = best_dofs
df["best_ngs_vta_volume"] = best_ngs_vta_volume
df["best_volume"] = best_volume
df.to_csv("vta_results_summary.csv", index=False)
