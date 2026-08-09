"""Aggregate benchmark records into a comparison table.

Reads every ``results/*.json`` written by ``run_benchmark.py`` and emits
``benchmark_summary.csv`` plus a Markdown table for the README.

Records are only comparable if they share ``benchmark_version`` and produce
the same problem size, so mismatches are reported rather than silently
averaged into the table.
"""

import argparse
import glob
import json
import os

import pandas as pd

RESULTS_DIR = "results"

# Timing columns, in pipeline order.
PHASES = [
    "load_images",
    "geometry",
    "contact_properties",
    "dielectric_model",
    "conductivity",
    "meshing_and_refinement_derived",
    "fem_solve",
    "point_model_copy",
    "time_reconstruction",
    "field_export",
]


def load_records(results_dir):
    """Load all benchmark records, newest first."""
    records = []
    for path in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(path) as fp:
            record = json.load(fp)
        record["_file"] = os.path.basename(path)
        records.append(record)
    return sorted(records, key=lambda r: r.get("timestamp_utc", ""), reverse=True)


def flatten(record):
    """Flatten one record into a single table row."""
    machine = record["machine"]
    best = record["best"]
    row = {
        "label": record["label"],
        "os": machine["os"],
        "arch": machine["arch"],
        "cpu_model": machine["cpu_model"],
        "cores_physical": machine["cores_physical"],
        "cores_logical": machine["cores_logical"],
        "memory_gb": machine["memory_gb"],
        "ngsolve": machine["packages"].get("ngsolve"),
        "neuron": machine["packages"].get("neuron"),
        "python": machine["python"],
        "commit": machine["ossdbs_commit"],
        "benchmark_version": record["benchmark_version"],
        "dofs": best["dofs"],
        "elements": best["elements"],
        "fem_total": best["fem_total"],
        "pam_total": best["pam_total"],
        "wall_total": best["wall_total"],
    }
    row.update({phase: best.get(phase) for phase in PHASES})
    return row


def check_comparability(df):
    """Warn about records that cannot be compared directly."""
    warnings = []
    for column, what in (
        ("benchmark_version", "benchmark version"),
        ("dofs", "problem size (DOFs)"),
    ):
        values = df[column].dropna().unique()
        if len(values) > 1:
            warnings.append(f"differing {what}: {sorted(values.tolist())}")
    return warnings


def to_markdown(df):
    """Render the headline comparison as a Markdown table.

    Written out by hand rather than via ``DataFrame.to_markdown``, which
    needs the optional ``tabulate`` package -- not worth adding a project
    dependency for a printed table.
    """
    columns = [
        "label",
        "os",
        "cpu_model",
        "cores_physical",
        "meshing_and_refinement_derived",
        "fem_solve",
        "fem_total",
        "pam_total",
        "wall_total",
    ]
    view = df[columns].rename(
        columns={
            "cores_physical": "cores",
            "meshing_and_refinement_derived": "mesh /s",
            "fem_solve": "solve /s",
            "fem_total": "FEM /s",
            "pam_total": "PAM /s",
            "wall_total": "total /s",
        }
    )
    headers = list(view.columns)
    rows = [["" if pd.isna(v) else str(v) for v in row] for row in view.to_numpy()]
    widths = [
        max([len(header), *(len(row[i]) for row in rows)])
        for i, header in enumerate(headers)
    ]

    def render(cells):
        padded = (c.ljust(w) for c, w in zip(cells, widths, strict=True))
        return "| " + " | ".join(padded) + " |"

    lines = [render(headers), "| " + " | ".join("-" * w for w in widths) + " |"]
    lines.extend(render(row) for row in rows)
    return "\n".join(lines)


def main():
    """Write the aggregated summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        default=RESULTS_DIR,
        help=f"Directory holding the benchmark JSON records (default: {RESULTS_DIR}).",
    )
    args = parser.parse_args()

    records = load_records(args.results_dir)
    if not records:
        print(f"No benchmark records found in {args.results_dir}/")
        return

    df = pd.DataFrame([flatten(record) for record in records])
    df.to_csv("benchmark_summary.csv", index=False)
    print(f"Wrote benchmark_summary.csv ({len(df)} record(s))")

    for warning in check_comparability(df):
        print(f"WARNING: {warning} -- these records are not directly comparable")

    print()
    print(to_markdown(df))


if __name__ == "__main__":
    main()
