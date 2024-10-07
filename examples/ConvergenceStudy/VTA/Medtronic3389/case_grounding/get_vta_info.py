import json
import os

import numpy as np
import pandas as pd

from ossdbs import VTAImage

vta_voxel_best = VTAImage("Results_VTA_best/VTA_solution_Lattice.nii")
best_volume = vta_voxel_best.get_vta_volume()
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
for result_dir in result_directories:
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm_info = json.load(fp)
    with open(os.path.join(result_dir, "Lattice.json")) as fp:
        lattice_info = json.load(fp)
    impedance_data = pd.read_csv(os.path.join(result_dir, "impedance.csv"))
    imp = impedance_data["real"].iloc[0]

    ngs_vta_volume = float(lattice_info["volume"])
    vta_voxel = VTAImage(os.path.join(result_dir, "VTA_solution_Lattice.nii"))
    dice = np.round(vta_voxel_best.compute_dice_coefficent(vta_voxel), 4)
    volume = np.round(vta_voxel.get_vta_volume(), 3)
    dofs = vcm_info["DOF"]
    time = vcm_info["Timings"]["ComputeSolution"][0]
    print(
        result_dir.replace("Results_VTA_", ""),
        " & ",
        dofs,
        " & ",
        np.round(time, 2),
        " & ",
        np.round(imp, 2),
        " & ",
        np.round(100.0 * abs(imp - best_imp) / best_imp, 2),
        " & ",
        np.round(ngs_vta_volume, 2),
        " & ",
        np.round(
            100.0 * abs(ngs_vta_volume - best_ngs_vta_volume) / best_ngs_vta_volume, 2
        ),
        " & ",
        volume,
        " & ",
        np.round(100.0 * abs(volume - best_volume) / best_volume, 2),
        " & ",
        dice,
        " \\\\",
    )
