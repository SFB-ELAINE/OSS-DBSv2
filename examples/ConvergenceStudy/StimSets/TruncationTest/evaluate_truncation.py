"""Collect per-protocol activations for the StimSets truncation study.

For every truncation ratio this writes ``<result_dir>_overview.csv`` with
one row per current protocol (in ``Current_protocols_0.csv`` order) and one
column per pathway -- the same layout produced by
``../evaluate_convergence.py`` for the mesh-convergence study.

The untruncated reference is by default the existing mesh-convergence run
``../Results_PAM_hp_material_refinement``, whose configuration is identical
to a truncation run with ``TruncateAfterActivePartRatio = None``. Its
overview CSV already exists, so nothing is recomputed for it. Pass
``--reference Results_PAM_truncation_none`` if you did re-run it here.
"""

import argparse
import json
import os

import pandas as pd

STIM_SETS_COMBINATIONS = 254

PATHWAYS = [
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

TRUNCATION_RATIOS = [5, 10, 20, 30]

# Untruncated run, reused from the mesh-convergence study by default.
DEFAULT_REFERENCE = os.path.join(os.path.pardir, "Results_PAM_hp_material_refinement")


def overview_path(result_dir):
    """Return the path of the overview CSV belonging to a result directory."""
    return result_dir + "_overview.csv"


def is_complete(result_dir):
    """Check whether a run wrote the status file of every pathway/protocol.

    The output directory is created when the FEM stage starts, long before
    the PAM sweep fills it, so directory existence alone is not enough --
    a still-running or interrupted run must be skipped, not half-read.
    """
    if not os.path.isdir(result_dir):
        return False
    last = STIM_SETS_COMBINATIONS - 1
    return all(
        os.path.isfile(
            os.path.join(result_dir, f"Pathway_status_{pathway}_{last}.json")
        )
        for pathway in PATHWAYS
    )


def build_overview(result_dir):
    """Read all Pathway_status JSONs of one run into a protocol x pathway table."""
    directory_results = {}
    for pathway in PATHWAYS:
        activations = []
        for stim_idx in range(STIM_SETS_COMBINATIONS):
            status_file = os.path.join(
                result_dir, f"Pathway_status_{pathway}_{stim_idx}.json"
            )
            with open(status_file) as fp:
                activations.append(json.load(fp)["percent_activated"])
        directory_results[pathway] = activations
    return pd.DataFrame(directory_results)


def main():
    """Write the overview CSVs for the reference and all truncated runs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        default=DEFAULT_REFERENCE,
        help=f"Untruncated reference result directory (default: {DEFAULT_REFERENCE}).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild overview CSVs even if they already exist.",
    )
    args = parser.parse_args()

    result_directories = [args.reference] + [
        f"Results_PAM_truncation_{ratio}" for ratio in TRUNCATION_RATIOS
    ]

    for result_dir in result_directories:
        out_file = overview_path(result_dir)
        if os.path.isfile(out_file) and not args.force:
            print(f"Skipping {out_file} (already exists, use --force to rebuild)")
            continue
        if not is_complete(result_dir):
            print(f"Skipping {result_dir} (missing or still running)")
            continue
        build_overview(result_dir).to_csv(out_file, index=False)
        print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
