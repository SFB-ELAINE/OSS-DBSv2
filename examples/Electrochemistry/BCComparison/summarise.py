"""Collect VTA volume, impedance, DOFs and timing from all six runs.

Mirrors examples/ConvergenceStudy/VTA/.../get_vta_info.py but reads
the BCComparison Results_* folders and writes a single CSV. Also
prints the per-contact floating potentials for the floating runs
(read from the per-run floating_potentials.csv, written by the FEM
when at least one contact is floating).
"""

import json
import os

import pandas as pd

PANEL_A = ["A1_floating", "A2_insulating", "A3_mild_Z"]
PANEL_B = ["B1_dirichlet", "B2_mild_Z", "B3_strong_Z"]


def collect_one(run_id: str, panel: str) -> dict:
    """Read one Results_<run_id>/ folder and return a row dict."""
    result_dir = f"Results_{run_id}"
    with open(os.path.join(result_dir, "VCM_report.json")) as fp:
        vcm = json.load(fp)
    impedance = pd.read_csv(os.path.join(result_dir, "impedance.csv"))

    return {
        "panel": panel,
        "run_id": run_id,
        "dofs": vcm["DOF"],
        "elements": vcm["Elements"],
        "time_s": vcm["Timings"]["ComputeSolution"][0],
        "ngs_vta_volume_mm3": float(vcm["VTA_volume_mm3"]),
        "impedance_real_Ohm": impedance["real"].iloc[0],
        "impedance_imag_Ohm": impedance["imag"].iloc[0],
    }


def collect_floating_potentials(run_ids: list[str]) -> pd.DataFrame:
    """Long-format table of floating contact voltages per run, single freq."""
    rows = []
    for run_id in run_ids:
        path = os.path.join(f"Results_{run_id}", "floating_potentials.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        # Single-frequency export: take the first row.
        row0 = df.iloc[0]
        contact_cols = [c for c in df.columns if c.endswith("_real")]
        for real_col in contact_cols:
            contact = real_col[: -len("_real")]
            imag_col = f"{contact}_imag"
            rows.append(
                {
                    "run_id": run_id,
                    "freq_Hz": float(row0["freq"]),
                    "contact": contact,
                    "V_real": float(row0[real_col]),
                    "V_imag": float(row0[imag_col]) if imag_col in df.columns else 0.0,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    """Aggregate all six runs into bc_comparison_summary.csv."""
    rows = [collect_one(rid, "A") for rid in PANEL_A]
    rows += [collect_one(rid, "B") for rid in PANEL_B]
    df = pd.DataFrame(rows)
    df.to_csv("bc_comparison_summary.csv", index=False)
    pd.set_option("display.float_format", lambda x: f"{x:.4g}")
    print(df.to_string(index=False))

    floating_df = collect_floating_potentials(PANEL_A + PANEL_B)
    if floating_df.empty:
        print(
            "\n(no floating_potentials.csv found in any Results_*/ — "
            "no floating contacts in any run)"
        )
    else:
        floating_df.to_csv("bc_comparison_floating_potentials.csv", index=False)
        print("\nFloating-contact voltages (single frequency):")
        print(floating_df.to_string(index=False))


if __name__ == "__main__":
    main()
