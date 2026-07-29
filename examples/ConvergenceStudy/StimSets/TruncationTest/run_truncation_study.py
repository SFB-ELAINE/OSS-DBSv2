# Truncation study: how many time points does NEURON need for a reliable
# activation estimate, evaluated over the full StimSets protocol sweep?
#
# This is the multi-protocol counterpart of
#   PAM_3/BostonVerciseDirected/TruncationTest/run_convergence_study.py
# All runs use the hp_material_refinement mesh strategy; the only thing
# varied is TruncateAfterActivePartRatio.
#
# In StimSets mode one run means eight per-contact unit solves
# (Results_PAM_truncation_<label>E1C1 ... E1C8) followed by a single
# run_PAM call that superimposes them for all 254 current protocols.
#
# Usage:
#   python run_truncation_study.py                 # run all ratios sequentially
#   python run_truncation_study.py 5 10            # run selected ratios
#   python run_truncation_study.py none            # run the untruncated reference
#   python run_truncation_study.py --list          # show ratios and disk estimates
#
# Each ratio runs in its own subprocess by default (disable with --no-isolated):
# StimSets plus HP refinement accumulates a lot of parent-process state, and a
# fresh interpreter per ratio keeps the peak RAM bounded. This is still strictly
# sequential -- one ratio finishes completely before the next starts.

import argparse
import json
import logging
import os
import subprocess
import sys
from copy import deepcopy

import ossdbs
from ossdbs.main import main_run

# Truncation ratios to be checked. ``None`` means no truncation and serves as
# the reference. The untruncated run is identical to the mesh-convergence
# strategy ``hp_material_refinement``, so the evaluation reuses that existing
# run by default -- see DEFAULT_REFERENCE in evaluate_truncation.py.
TRUNCATION_RATIOS = [None, 5, 10, 20, 30]

# Study directory sits one level below the StimSets root, which holds
# oss-dbs_parameters.json, segmask.nii.gz and OSS_sim_files_rh/.
STUDY_ROOT = os.path.pardir


def ratio_label(ratio):
    """Return the directory label used for a truncation ratio."""
    return "none" if ratio is None else str(ratio)


def output_path(ratio):
    """Return the OutputPath for a truncation ratio."""
    return f"Results_PAM_truncation_{ratio_label(ratio)}"


def remove_file_handler(logger):
    """Remove file handler from logger instance."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


def save_input_dict(base_input_dict):
    """Save input dictionary for the truncation run.

    The saved dict lives inside the output directory, i.e. one level
    deeper than this script, so every input path gains one more ``..``.
    """
    input_dict = deepcopy(base_input_dict)
    for key_path in (
        ("MaterialDistribution", "MRIPath"),
        ("PointModel", "Pathway", "FileName"),
        ("StimSets", "StimSetsFile"),
    ):
        node = input_dict
        for key in key_path[:-1]:
            node = node[key]
        node[key_path[-1]] = os.path.join(os.path.pardir, node[key_path[-1]])
    input_dict["PathwayFile"] = os.path.join(os.path.pardir, input_dict["PathwayFile"])
    input_dict["OutputPath"] = "./"

    with open(
        os.path.join(base_input_dict["OutputPath"], "input_dict.json"), "w"
    ) as fp:
        json.dump(input_dict, fp, indent=2)


def run_strategy(input_dict, logger, stage="all"):
    """Run the FEM unit solves and/or the StimSets PAM sweep.

    ``stage="fem"`` stops after the eight unit solutions, ``stage="pam"``
    reuses unit solutions already on disk. Splitting the two is useful to
    validate the unit solutions before committing to the multi-hour PAM
    sweep, and to resume a run whose PAM step was interrupted.
    """
    if stage in ("all", "fem"):
        main_run(input_dict)
        save_input_dict(input_dict)
    if stage in ("all", "pam"):
        ossdbs.api.run_PAM(input_dict)
    remove_file_handler(logger)


def setup_base_config():
    """Load the StimSets configuration and clean its mesh settings."""
    with open(os.path.join(STUDY_ROOT, "oss-dbs_parameters.json")) as fp:
        base = json.load(fp)

    # adjust paths: inputs live in the StimSets root, one level up
    base["MaterialDistribution"]["MRIPath"] = os.path.join(
        STUDY_ROOT, base["MaterialDistribution"]["MRIPath"]
    )
    base["PointModel"]["Pathway"]["FileName"] = os.path.join(
        STUDY_ROOT, base["PointModel"]["Pathway"]["FileName"]
    )
    base["PathwayFile"] = os.path.join(STUDY_ROOT, base["PathwayFile"])
    base["StimSets"]["StimSetsFile"] = os.path.join(
        STUDY_ROOT, base["StimSets"]["StimSetsFile"]
    )

    # do not scale solution, the current protocols carry the amplitudes
    base["Scaling"] = 1.0

    # clean base state
    base["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    base["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
    base["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    base["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = ""
    base["Mesh"]["MaterialRefinementSteps"] = 0

    return base


def configure_hp_material_refinement(base):
    """Default mesh + HP ref. + 1x material ref."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {
        "Active": True,
        "Levels": 2,
        "Factor": 0.125,
    }
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    return cfg


def build_config(ratio):
    """Return the full configuration for one truncation ratio."""
    cfg = configure_hp_material_refinement(setup_base_config())
    if ratio is not None:
        cfg["TruncateAfterActivePartRatio"] = float(ratio)
    cfg["OutputPath"] = output_path(ratio)
    return cfg


def run_single(ratio, loglevel, stage="all"):
    """Run one truncation ratio in this process."""
    ossdbs.set_logger(level=loglevel)
    logger = logging.getLogger("ossdbs")
    run_strategy(build_config(ratio), logger, stage=stage)


def run_isolated(ratio, loglevel, stage="all"):
    """Run one truncation ratio in a fresh interpreter."""
    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--single",
        ratio_label(ratio),
        "--stage",
        stage,
        "--loglevel",
        str(loglevel),
    ]
    subprocess.run(cmd, check=True, cwd=os.path.dirname(os.path.abspath(__file__)))


def parse_ratio(token):
    """Resolve a CLI token to a truncation ratio."""
    if token.lower() == "none":
        return None
    try:
        value = int(token)
    except ValueError:
        return "invalid"
    return value if value in TRUNCATION_RATIOS else "invalid"


def resolve_ratios(tokens, parser):
    """Resolve CLI tokens to truncation ratios, or return all of them."""
    if not tokens:
        return list(TRUNCATION_RATIOS)
    ratios = []
    for token in tokens:
        ratio = parse_ratio(token)
        if ratio == "invalid":
            parser.error(
                f"Unknown truncation ratio: {token!r}. Use --list to see options."
            )
        ratios.append(ratio)
    return ratios


def print_ratio_table():
    """Print the truncation ratios with the kept fraction of the signal."""
    signal = setup_base_config()["StimulationSignal"]
    active_us = (
        signal["PulseWidth[us]"]
        + signal["InterPulseWidth[us]"]
        + signal["CounterPulseWidth[us]"]
    )
    period_us = 1e6 / signal["Frequency[Hz]"]
    print(f"Active signal part: {active_us:g} us, period {period_us:.0f} us")
    print(f"{'ratio':>6s}  {'kept signal':>12s}  {'fraction':>9s}")
    for ratio in TRUNCATION_RATIOS:
        if ratio is None:
            print(f"{'none':>6s}  {period_us:11.0f} us  {1.0:8.1%}")
        else:
            kept = ratio * active_us
            print(f"{ratio:>6d}  {kept:11.0f} us  {kept / period_us:8.1%}")


def main():
    """Run the selected truncation ratios sequentially."""
    parser = argparse.ArgumentParser(
        description="Run the StimSets PAM truncation study.",
    )
    parser.add_argument(
        "ratios",
        nargs="*",
        help="Truncation ratios to run (5, 10, 20, 30 or none). Default: all.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List the ratios with signal length and disk estimates, then exit.",
    )
    parser.add_argument(
        "--single",
        metavar="RATIO",
        help="Internal: run exactly this ratio in the current process.",
    )
    parser.add_argument(
        "--stage",
        choices=("all", "fem", "pam"),
        default="all",
        help=(
            "Which part to run: 'fem' only the eight unit solves, "
            "'pam' only the protocol sweep on existing unit solves "
            "(default: all)."
        ),
    )
    parser.add_argument(
        "--no-isolated",
        dest="isolated",
        action="store_false",
        help="Run all ratios in this process instead of one subprocess each.",
    )
    parser.add_argument(
        "--loglevel",
        type=int,
        default=logging.INFO,
        help="Logging verbosity (10=DEBUG, 20=INFO).",
    )
    args = parser.parse_args()

    if args.list:
        print_ratio_table()
        return

    if args.single is not None:
        ratio = parse_ratio(args.single)
        if ratio == "invalid":
            parser.error(f"Unknown truncation ratio: {args.single!r}")
        run_single(ratio, args.loglevel, stage=args.stage)
        return

    to_run = resolve_ratios(args.ratios, parser)

    for ratio in to_run:
        label = ratio_label(ratio)
        print(f"\n{'=' * 60}")
        print(f"Running truncation ratio: {label} (stage: {args.stage})")
        print(f"{'=' * 60}")
        if args.isolated:
            run_isolated(ratio, args.loglevel, stage=args.stage)
        else:
            run_single(ratio, args.loglevel, stage=args.stage)

    print(f"\nCompleted {len(to_run)} truncation ratio(s).")


if __name__ == "__main__":
    main()
