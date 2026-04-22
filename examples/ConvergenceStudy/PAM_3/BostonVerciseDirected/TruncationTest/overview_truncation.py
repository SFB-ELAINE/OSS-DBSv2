"""Generate truncation_results_summary.csv for the truncation test.

Compares per-pathway activation at different TruncateAfterActivePartRatio
values against the untruncated run (Results_PAM_truncation_none).
"""

import json
import os

import h5py
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

ref_dir = "Results_PAM_truncation_none"
truncation_ratios = [5, 10, 20, 30]

# Load reference activations
ref_activations = {}
for pw in pathways:
    with open(os.path.join(ref_dir, f"Pathway_status_{pw}.json")) as f:
        ref_activations[pw] = json.load(f)["percent_activated"]

rows = []
# All directories to process: reference first, then truncated
all_dirs = [(ref_dir, "none")] + [
    (f"Results_PAM_truncation_{r}", str(r)) for r in truncation_ratios
]

for result_dir, label in all_dirs:
    # Time steps from H5
    h5_path = os.path.join(result_dir, "oss_time_result_PAM.h5")
    with h5py.File(h5_path, "r") as f:
        timesteps = f[f"{pathways[0]}/axon0/Potential[V]"].shape[1]

    # H5 file size in MB
    h5_size_mb = os.path.getsize(h5_path) / 1e6

    # VCM timings
    with open(os.path.join(result_dir, "VCM_report.json")) as f:
        vcm = json.load(f)
    dofs = vcm["DOF"]
    vcm_time = sum(vcm["Timings"]["ComputeSolution"])

    row = {
        "truncation_ratio": label,
        "timesteps": timesteps,
        "h5_size_MB": round(h5_size_mb, 1),
        "dofs": dofs,
        "vcm_time": round(vcm_time, 1),
    }

    for pw in pathways:
        pw_file = os.path.join(result_dir, f"Pathway_status_{pw}.json")
        with open(pw_file) as f:
            activated = json.load(f)["percent_activated"]
        error = abs(activated - ref_activations[pw])
        row[f"{pw}_activated"] = activated
        row[f"{pw}_error"] = error

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("truncation_results_summary.csv", index=False)
print("Wrote truncation_results_summary.csv")
print(df[["truncation_ratio", "timesteps", "h5_size_MB"]].to_string(index=False))
print()
for pw in pathways:
    vals = df[["truncation_ratio", f"{pw}_activated", f"{pw}_error"]].to_string(
        index=False
    )
    print(f"{pw}:")
    print(vals)
    print()
