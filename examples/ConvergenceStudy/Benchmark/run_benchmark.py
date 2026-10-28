"""Benchmark the best convergence-study strategy on the current machine.

Strategy under test (pinned here, not imported from a study script, so the
benchmark stays comparable as the convergence studies evolve):

    Fine mesh + HP refinement (2 levels, factor 0.125) + 1x material refinement

Case: PAM_3 / BostonVerciseDirected, single monopolar current-controlled
protocol. Small enough to run on a laptop, which is the point -- the
benchmark is only useful if it is run on more than one machine.

Both stages are timed and both are required, so the benchmark needs NEURON
and therefore runs on Linux and macOS only.

Usage:
    python run_benchmark.py                 # one run, writes results/<...>.json
    python run_benchmark.py --repeats 3     # three runs, keeps the fastest
    python run_benchmark.py --label "hpc-node-01"
    python run_benchmark.py --dry-run       # print the machine context and exit
"""

import argparse
import json
import logging
import os
import platform
import shutil
import socket
import time
from copy import deepcopy
from datetime import datetime, timezone

import machine_info

import ossdbs
from ossdbs.main import main_run

# Inputs live in the sibling PAM_3/ directory.
STUDY_ROOT = os.path.join(os.path.pardir, "PAM_3")
# Benchmark/ -> ConvergenceStudy/ -> examples/ -> repo root
REPO_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), *([os.path.pardir] * 3)
)
RESULTS_DIR = "results"
OUTPUT_PATH = "Results_benchmark"

# The pinned strategy. Changing any of this invalidates comparability with
# previously recorded results, so bump BENCHMARK_VERSION when you do.
BENCHMARK_VERSION = 1
MESHING_HYPOTHESIS = "Fine"
HP_REFINEMENT = {"Active": True, "Levels": 2, "Factor": 0.125}
MATERIAL_REFINEMENT_STEPS = 1


def build_config():
    """Return the pinned benchmark configuration."""
    with open(os.path.join(STUDY_ROOT, "oss-dbs_parameters.json")) as fp:
        cfg = json.load(fp)

    for key_path in (
        ("MaterialDistribution", "MRIPath"),
        ("PointModel", "Pathway", "FileName"),
    ):
        node = cfg
        for key in key_path[:-1]:
            node = node[key]
        node[key_path[-1]] = os.path.join(STUDY_ROOT, node[key_path[-1]])
    cfg["PathwayFile"] = os.path.join(STUDY_ROOT, cfg["PathwayFile"])

    # single protocol, unscaled
    cfg["Scaling"] = 1.0
    cfg["ScalingIndex"] = None

    # pinned mesh strategy
    cfg["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = MESHING_HYPOTHESIS
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    cfg["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = ""
    cfg["Mesh"]["HPRefinement"] = deepcopy(HP_REFINEMENT)
    cfg["Mesh"]["MaterialRefinementSteps"] = MATERIAL_REFINEMENT_STEPS

    cfg["OutputPath"] = OUTPUT_PATH
    return cfg


def _total(value):
    """Sum a timing entry that may be per-frequency or a single number.

    The volume conductor records some phases once per frequency
    (``ComputeSolution``, ``CopyValues``) and others once per run
    (``FieldExport``), so both shapes have to be handled.
    """
    if value is None:
        return 0.0
    if isinstance(value, list):
        return float(sum(value))
    return float(value)


def phase_timings():
    """Read the persisted phase timings of the last run.

    ``run_report.json`` holds the top-level phases recorded by ``main_run``;
    ``VCM_report.json`` holds the FEM-internal breakdown plus problem size.
    """
    with open(os.path.join(OUTPUT_PATH, "run_report.json")) as fp:
        run_report = json.load(fp)["Timings"]
    with open(os.path.join(OUTPUT_PATH, "VCM_report.json")) as fp:
        vcm_report = json.load(fp)

    vcm = vcm_report["Timings"]
    solve = _total(vcm.get("ComputeSolution"))
    copy_values = _total(vcm.get("CopyValues"))
    reconstruct = sum(
        _total(value)
        for key, value in vcm.items()
        if key.startswith("ReconstructTimeSignals_PointModel")
    )
    field_export = _total(vcm.get("FieldExport"))

    # main_run's "VolumeConductor" bucket covers mesh generation, refinement,
    # FE space setup and the full analysis. Subtracting the analysis timings
    # isolates meshing + refinement, which is the part that scales with
    # single-core speed rather than with core count -- worth separating when
    # comparing machines.
    analysis = solve + copy_values + reconstruct + field_export
    meshing = run_report.get("VolumeConductor", 0.0) - analysis

    return {
        "load_images": round(run_report.get("MRI", 0.0), 3),
        "geometry": round(run_report.get("ModelGeometry", 0.0), 3),
        "contact_properties": round(run_report.get("ContactProperties", 0.0), 3),
        "dielectric_model": round(run_report.get("DielectricModel", 0.0), 3),
        "conductivity": round(run_report.get("ConductivityCF", 0.0), 3),
        "meshing_and_refinement_derived": round(meshing, 3),
        "fem_solve": round(solve, 3),
        "point_model_copy": round(copy_values, 3),
        "time_reconstruction": round(reconstruct, 3),
        "field_export": round(field_export, 3),
        "n_frequencies": len(vcm["ComputeSolution"]),
        "dofs": vcm_report["DOF"],
        "elements": vcm_report["Elements"],
    }


def remove_file_handler(logger):
    """Remove file handler so repeated runs do not stack log handlers."""
    for handler in list(logger.handlers):
        if isinstance(handler, logging.FileHandler):
            logger.removeHandler(handler)


def single_run(loglevel):
    """Run FEM + PAM once and return the timing record."""
    if os.path.isdir(OUTPUT_PATH):
        shutil.rmtree(OUTPUT_PATH)

    ossdbs.set_logger(level=loglevel)
    logger = logging.getLogger("ossdbs")
    cfg = build_config()

    fem_start = time.perf_counter()
    main_run(cfg)
    fem_total = time.perf_counter() - fem_start

    pam_start = time.perf_counter()
    ossdbs.api.run_PAM(cfg)
    pam_total = time.perf_counter() - pam_start

    remove_file_handler(logger)

    record = {
        "fem_total": round(fem_total, 3),
        "pam_total": round(pam_total, 3),
        "wall_total": round(fem_total + pam_total, 3),
    }
    record.update(phase_timings())
    return record


def main():
    """Run the benchmark and record the result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Number of runs; the fastest by wall_total is recorded (default: 1).",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Name for this machine in the results table (default: hostname).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the machine context and exit without running.",
    )
    parser.add_argument(
        "--loglevel",
        type=int,
        default=logging.INFO,
        help="Logging verbosity (10=DEBUG, 20=INFO).",
    )
    args = parser.parse_args()

    context = machine_info.collect(REPO_ROOT)

    if args.dry_run:
        print(json.dumps(context, indent=2))
        return

    if platform.system() == "Windows":
        parser.error(
            "The benchmark times FEM and PAM together and PAM needs NEURON, "
            "which is not available on Windows."
        )
    if context["packages"]["neuron"] is None:
        parser.error("NEURON is not importable; install it to run the PAM stage.")

    runs = []
    for index in range(args.repeats):
        print(f"\n{'=' * 60}\nBenchmark run {index + 1}/{args.repeats}\n{'=' * 60}")
        runs.append(single_run(args.loglevel))
        print(f"wall_total: {runs[-1]['wall_total']} s")

    best = min(runs, key=lambda run: run["wall_total"])
    label = args.label or socket.gethostname()
    result = {
        "benchmark_version": BENCHMARK_VERSION,
        "label": label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "strategy": {
            "meshing_hypothesis": MESHING_HYPOTHESIS,
            "hp_refinement": HP_REFINEMENT,
            "material_refinement_steps": MATERIAL_REFINEMENT_STEPS,
            "case": "PAM_3/BostonVerciseDirected, monopolar single protocol",
        },
        "machine": context,
        "repeats": args.repeats,
        "best": best,
        "all_runs": runs,
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    stamp = result["timestamp_utc"].replace(":", "").replace("-", "")
    safe_label = "".join(c if c.isalnum() or c in "-_" else "-" for c in label)
    out_file = os.path.join(RESULTS_DIR, f"{safe_label}-{stamp}.json")
    with open(out_file, "w") as fp:
        json.dump(result, fp, indent=2)

    print(f"\nWrote {out_file}")
    print(f"  DOFs {best['dofs']}, elements {best['elements']}")
    print(f"  FEM  {best['fem_total']:8.1f} s")
    print(f"  PAM  {best['pam_total']:8.1f} s")
    print(f"  total{best['wall_total']:8.1f} s")


if __name__ == "__main__":
    main()
