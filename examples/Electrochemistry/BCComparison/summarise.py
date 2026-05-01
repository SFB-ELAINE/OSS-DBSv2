"""Collect VTA volume, impedance, DOFs and timing from all six runs.

Mirrors examples/ConvergenceStudy/VTA/.../get_vta_info.py but reads
the BCComparison Results_* folders and writes a single CSV.
"""

import json
import os

import pandas as pd

from ossdbs import VTAImage

PANEL_A = ["A1_floating", "A2_insulating", "A3_mild_Z"]
PANEL_B = ["B1_dirichlet", "B2_mild_Z", "B3_strong_Z"]


def collect_one(run_id: str, panel: str) -> dict:
    """Read one Results_<run_id>/ folder and return a row dict."""
    result_dir = f"Results_{run_id}"
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm = json.load(fp)
    with open(os.path.join(result_dir, "Lattice.json")) as fp:
        lattice = json.load(fp)
    impedance = pd.read_csv(os.path.join(result_dir, "impedance.csv"))
    voxel = VTAImage(os.path.join(result_dir, "VTA_solution_Lattice.nii.gz"))

    return {
        "panel": panel,
        "run_id": run_id,
        "dofs": vcm["DOF"],
        "elements": vcm["Elements"],
        "time_s": vcm["Timings"]["ComputeSolution"][0],
        "ngs_vta_volume_mm3": float(lattice["volume"]),
        "voxel_vta_volume_mm3": voxel.get_vta_volume(),
        "impedance_real_Ohm": impedance["real"].iloc[0],
        "impedance_imag_Ohm": impedance["imag"].iloc[0],
    }


def main() -> None:
    """Aggregate all six runs into bc_comparison_summary.csv."""
    rows = [collect_one(rid, "A") for rid in PANEL_A]
    rows += [collect_one(rid, "B") for rid in PANEL_B]
    df = pd.DataFrame(rows)
    df.to_csv("bc_comparison_summary.csv", index=False)
    pd.set_option("display.float_format", lambda x: f"{x:.4g}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
