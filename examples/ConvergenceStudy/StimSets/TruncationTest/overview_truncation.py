"""Generate truncation_results_summary.csv for the StimSets truncation study.

Multi-protocol counterpart of
``PAM_3/BostonVerciseDirected/TruncationTest/overview_truncation.py``.

Where the single-protocol study compares one activation value per pathway,
this compares the whole StimSets sweep: for every truncation ratio and every
pathway it reports the mean activation over the current protocols, plus the
mean and maximum absolute deviation from the untruncated reference across
those protocols. The maximum is the interesting number here -- a truncation
ratio is only safe if the *worst* protocol still agrees.

Only clinically feasible protocols are aggregated: total current
sum(|I|) <= L1_THRESHOLD_MA and every contact within
[LOWER_LIM_MA, UPPER_LIM_MA], matching ``../overview_refinement.py``.
The overview CSVs stay full (one row per protocol, in
``Current_protocols_0.csv`` order); the filter selects rows positionally.

Cost columns are per-run totals over the eight per-contact unit solves
(``...E1C1`` ... ``...E1C8``), since a StimSets run needs all of them.

Run ``evaluate_truncation.py`` first to produce the overview CSVs.
"""

import argparse
import json
import os

import h5py
import numpy as np
import pandas as pd
from evaluate_truncation import (
    DEFAULT_REFERENCE,
    PATHWAYS,
    TRUNCATION_RATIOS,
    overview_path,
)

# Clinical-feasibility filter on the current protocols (see module docstring).
STIM_SETS_FILE = os.path.join(
    os.path.pardir, "OSS_sim_files_rh", "Current_protocols_0.csv"
)
L1_THRESHOLD_MA = 8.0
LOWER_LIM_MA = -4.0
UPPER_LIM_MA = 4.0

N_CONTACTS = 8


def kept_protocol_indices():
    """Positional indices of protocols passing the clinical-feasibility filter.

    The overview CSVs share the row order of ``Current_protocols_0.csv``, so
    these indices select the same protocols there via ``.iloc``.
    """
    arr = np.nan_to_num(pd.read_csv(STIM_SETS_FILE).to_numpy(dtype=float), nan=0.0)
    within_l1 = np.abs(arr).sum(axis=1) <= L1_THRESHOLD_MA
    within_limits = ((arr >= LOWER_LIM_MA) & (arr <= UPPER_LIM_MA)).all(axis=1)
    return np.flatnonzero(within_l1 & within_limits)


def contact_dirs(result_dir):
    """Return the existing per-contact unit-solution directories of a run."""
    dirs = [f"{result_dir}E1C{i}" for i in range(1, N_CONTACTS + 1)]
    return [d for d in dirs if os.path.isdir(d)]


def run_costs(result_dir):
    """Return time steps, H5 size, DOFs and solve time of one StimSets run."""
    dirs = contact_dirs(result_dir)
    if not dirs:
        raise FileNotFoundError(f"No per-contact directories found for {result_dir}")

    h5_paths = [os.path.join(d, "oss_time_result_PAM.h5") for d in dirs]
    with h5py.File(h5_paths[0], "r") as f:
        timesteps = f[f"{PATHWAYS[0]}/axon0/Potential[V]"].shape[1]

    with open(os.path.join(dirs[0], "VCM_report.json")) as f:
        dofs = json.load(f)["DOF"]

    vcm_time = 0.0
    for d in dirs:
        with open(os.path.join(d, "VCM_report.json")) as f:
            vcm_time += sum(json.load(f)["Timings"]["ComputeSolution"])

    return {
        "timesteps": timesteps,
        "h5_size_MB": round(os.path.getsize(h5_paths[0]) / 1e6, 1),
        "h5_size_total_MB": round(sum(os.path.getsize(p) for p in h5_paths) / 1e6, 1),
        "n_contacts": len(dirs),
        "dofs": dofs,
        "vcm_time": round(vcm_time, 1),
    }


def main():
    """Write truncation_results_summary.csv."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        default=DEFAULT_REFERENCE,
        help=f"Untruncated reference result directory (default: {DEFAULT_REFERENCE}).",
    )
    args = parser.parse_args()

    kept_idx = kept_protocol_indices()
    print(f"Aggregating over {len(kept_idx)} clinically feasible protocols")

    ref_df = (
        pd.read_csv(overview_path(args.reference)).iloc[kept_idx].reset_index(drop=True)
    )

    all_runs = [(args.reference, "none")] + [
        (f"Results_PAM_truncation_{ratio}", str(ratio)) for ratio in TRUNCATION_RATIOS
    ]

    rows = []
    for result_dir, label in all_runs:
        if not os.path.isfile(overview_path(result_dir)):
            print(f"Skipping ratio {label} ({overview_path(result_dir)} missing)")
            continue

        row = {"truncation_ratio": label}
        row.update(run_costs(result_dir))

        df = (
            pd.read_csv(overview_path(result_dir)).iloc[kept_idx].reset_index(drop=True)
        )
        for pathway in PATHWAYS:
            abs_errors = np.abs(df[pathway].to_numpy() - ref_df[pathway].to_numpy())
            row[f"{pathway}_activated"] = df[pathway].mean()
            row[f"{pathway}_error"] = abs_errors.mean()
            row[f"{pathway}_max_error"] = abs_errors.max()
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv("truncation_results_summary.csv", index=False)
    print("Wrote truncation_results_summary.csv")
    print(
        summary[
            [
                "truncation_ratio",
                "timesteps",
                "h5_size_MB",
                "h5_size_total_MB",
                "dofs",
                "vcm_time",
            ]
        ].to_string(index=False)
    )
    print()
    for pathway in PATHWAYS:
        cols = [
            "truncation_ratio",
            f"{pathway}_activated",
            f"{pathway}_error",
            f"{pathway}_max_error",
        ]
        print(f"{pathway}:")
        print(summary[cols].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
