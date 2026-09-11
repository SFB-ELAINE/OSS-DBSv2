"""Benchmark the best convergence-study strategy on a VTA workload.

Companion to ``run_benchmark.py``. The pinned mesh strategy is identical,
so the two records isolate the effect of the workload rather than of the
refinement:

    Fine mesh + HP refinement (2 levels, factor 0.125) + 1x material refinement

Case: VTA / BostonScientificVerciseDirected, case grounding -- contact 2 at
1 V against the brain surface at 0 V. Same electrode as the PAM benchmark,
monopolar, which makes it the direct analogue. The other configurations in
``../VTA/BostonScientificVerciseDirected`` (adjacent and distant contacts)
are bipolar variants.

The VTA workload differs from the PAM one in ways that show up in the
breakdown:

* voltage-controlled, single-frequency Multisine -- one solve instead of
  twelve, and no time-domain reconstruction at all
* an 864000-point lattice instead of 75935 pathway points, so locating
  points on the mesh dominates rather than the solve
* no NEURON, so unlike the PAM benchmark this one also runs on Windows

Usage:
    python run_benchmark_vta.py                 # one run
    python run_benchmark_vta.py --repeats 3     # keep the fastest
    python run_benchmark_vta.py --label hpc-node-01
    python run_benchmark_vta.py --dry-run       # machine context only
"""

import argparse
import json
import logging
import os
import shutil
import socket
import time
from copy import deepcopy
from datetime import datetime, timezone

import machine_info

import ossdbs
from ossdbs.main import main_run

STUDY_ROOT = os.path.join(os.path.pardir, "VTA")
# the VTA segmentation mask is byte-identical to the PAM_3 one, which is the
# copy tracked in the repository, so point at that instead of requiring a
# second 1.1 MB binary
MRI_PATH = os.path.join(os.path.pardir, "PAM_3", "segmask.nii.gz")
# Benchmark/ -> ConvergenceStudy/ -> examples/ -> repo root
REPO_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), *([os.path.pardir] * 3)
)
# Two variants of the same case, differing only in how the VTA is obtained.
#
# "lattice"  samples the solution on an 864000-point grid, writes the CSV and
#            the E-field / VTA NIfTI images, and additionally integrates the
#            VTA volume on the FEM mesh.
# "ngsolve"  drops the lattice entirely and takes the VTA volume purely from
#            the NGSolve integration over the mesh. No point evaluation, no
#            NIfTI export.
#
# Both produce the same vta_volume_mm3, because that number never came from
# the lattice in the first place -- it is computed by
# threshold_frequency_domain_Efield on the FEM mesh. Comparing the two
# therefore prices the lattice sampling and image export, which is the bulk
# of the "lattice" runtime, without changing the physics being measured.
VARIANTS = {
    "lattice": {
        "results_dir": "results_vta",
        "output_path": "Results_benchmark_vta",
        "case": "VTA/BostonScientificVerciseDirected, case grounding",
    },
    "ngsolve": {
        "results_dir": "results_vta_ngsolve",
        "output_path": "Results_benchmark_vta_ngsolve",
        "case": (
            "VTA/BostonScientificVerciseDirected, case grounding "
            "(NGSolve VTA only, no lattice or NIfTI export)"
        ),
    },
}

# Pinned strategy and case. Raise BENCHMARK_VERSION whenever either changes,
# so collect_results.py flags older records instead of mixing them in.
BENCHMARK_VERSION = 1
MESHING_HYPOTHESIS = "Fine"
HP_REFINEMENT = {"Active": True, "Levels": 2, "Factor": 0.125}
MATERIAL_REFINEMENT_STEPS = 1
ELECTRODE_NAME = "BostonScientificVerciseDirected"
ACTIVE_CONTACT_ID = 2
ACTIVE_CONTACT_VOLTAGE = 1.0

BASE_CONTACT = {
    "Contact_ID": 1,
    "Active": False,
    "Current[A]": 0.0,
    "Voltage[V]": 0.0,
    "Floating": False,
}


def build_config(variant="lattice"):
    """Return the pinned VTA benchmark configuration for one variant."""
    with open(os.path.join(STUDY_ROOT, "base_settings.json")) as fp:
        cfg = json.load(fp)

    cfg["MaterialDistribution"]["MRIPath"] = MRI_PATH

    # electrode with a single active contact
    cfg["Electrodes"][0]["Name"] = ELECTRODE_NAME
    contact = deepcopy(BASE_CONTACT)
    contact["Contact_ID"] = ACTIVE_CONTACT_ID
    contact["Active"] = True
    contact["Voltage[V]"] = ACTIVE_CONTACT_VOLTAGE
    cfg["Electrodes"][0]["Contacts"].append(contact)

    # case grounding: brain surface held at 0 V
    cfg["Surfaces"] = [{"Name": "BrainSurface", "Active": True, "Voltage[V]": 0.0}]

    # lattice of the case_grounding study; switched off entirely in the
    # ngsolve variant, which removes the point evaluation and with it the
    # CSV and NIfTI exports that are driven by the point model
    cfg["PointModel"]["Lattice"]["Shape"] = {"x": 60, "y": 60, "z": 240}
    cfg["PointModel"]["Lattice"]["PointDistance[mm]"] = 0.125
    cfg["PointModel"]["Lattice"]["Active"] = variant == "lattice"

    # pinned mesh strategy, matching run_benchmark.py
    cfg["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = MESHING_HYPOTHESIS
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    cfg["Mesh"]["HPRefinement"] = deepcopy(HP_REFINEMENT)
    cfg["Mesh"]["MaterialRefinementSteps"] = MATERIAL_REFINEMENT_STEPS

    cfg["OutputPath"] = VARIANTS[variant]["output_path"]
    return cfg


def _total(value):
    """Sum a timing entry that may be per-frequency or a single number."""
    if value is None:
        return 0.0
    if isinstance(value, list):
        return float(sum(value))
    return float(value)


def phase_timings(output_path):
    """Read the persisted phase timings of the last run."""
    with open(os.path.join(output_path, "run_report.json")) as fp:
        run_report = json.load(fp)["Timings"]
    with open(os.path.join(output_path, "VCM_report.json")) as fp:
        vcm_report = json.load(fp)

    vcm = vcm_report["Timings"]
    solve = _total(vcm.get("ComputeSolution"))
    copy_values = _total(vcm.get("CopyValues"))
    field_export = _total(vcm.get("FieldExport"))
    # a single-frequency Multisine skips the time-domain reconstruction
    reconstruct = sum(
        _total(value)
        for key, value in vcm.items()
        if key.startswith("ReconstructTimeSignals_PointModel")
    )

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
        "n_frequencies": len(vcm.get("ComputeSolution", [])),
        "dofs": vcm_report["DOF"],
        "elements": vcm_report["Elements"],
        # not a timing, but the check that two machines solved the same
        # problem rather than merely spending the same time
        "vta_volume_mm3": vcm_report.get("VTA_volume_mm3"),
    }


def remove_file_handler(logger):
    """Remove file handler so repeated runs do not stack log handlers."""
    for handler in list(logger.handlers):
        if isinstance(handler, logging.FileHandler):
            logger.removeHandler(handler)


def single_run(loglevel, variant):
    """Run the VTA analysis once and return the timing record."""
    output_path = VARIANTS[variant]["output_path"]
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    ossdbs.set_logger(level=loglevel)
    logger = logging.getLogger("ossdbs")
    cfg = build_config(variant)

    start = time.perf_counter()
    main_run(cfg)
    fem_total = time.perf_counter() - start

    remove_file_handler(logger)

    record = {
        "fem_total": round(fem_total, 3),
        "pam_total": None,  # no NEURON stage in the VTA workload
        "wall_total": round(fem_total, 3),
    }
    record.update(phase_timings(output_path))
    return record


def main():
    """Run the VTA benchmark and record the result."""
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
        "--variant",
        choices=sorted(VARIANTS),
        default="lattice",
        help=(
            "'lattice' samples the 864000-point grid and writes the CSV and "
            "NIfTI exports; 'ngsolve' drops the lattice and takes the VTA "
            "volume from the NGSolve mesh integration alone (default: lattice)."
        ),
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

    runs = []
    for index in range(args.repeats):
        print(
            f"\n{'=' * 60}\nVTA benchmark ({args.variant}) "
            f"run {index + 1}/{args.repeats}\n{'=' * 60}"
        )
        runs.append(single_run(args.loglevel, args.variant))
        print(f"wall_total: {runs[-1]['wall_total']} s")

    best = min(runs, key=lambda run: run["wall_total"])
    label = args.label or socket.gethostname()
    result = {
        "benchmark_version": BENCHMARK_VERSION,
        "label": label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "variant": args.variant,
        "strategy": {
            "meshing_hypothesis": MESHING_HYPOTHESIS,
            "hp_refinement": HP_REFINEMENT,
            "material_refinement_steps": MATERIAL_REFINEMENT_STEPS,
            "case": VARIANTS[args.variant]["case"],
        },
        "machine": context,
        "repeats": args.repeats,
        "best": best,
        "all_runs": runs,
    }

    results_dir = VARIANTS[args.variant]["results_dir"]
    os.makedirs(results_dir, exist_ok=True)
    stamp = result["timestamp_utc"].replace(":", "").replace("-", "")
    safe_label = "".join(c if c.isalnum() or c in "-_" else "-" for c in label)
    out_file = os.path.join(results_dir, f"{safe_label}-{stamp}.json")
    with open(out_file, "w") as fp:
        json.dump(result, fp, indent=2)

    print(f"\nWrote {out_file}")
    print(f"  DOFs {best['dofs']}, elements {best['elements']}")
    print(f"  VTA volume {best['vta_volume_mm3']} mm^3")
    print(f"  total{best['wall_total']:8.1f} s")


if __name__ == "__main__":
    main()
