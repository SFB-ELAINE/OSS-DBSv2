# Convergence study: mesh refinement strategies for StimSets VTA
#
# Same StimSets unit-solution workflow as the PAM study, but with a
# Lattice point model instead of a Pathway. Each strategy runs the
# per-contact unit solves (run_stim_sets) and writes one
# `E_field_Lattice.csv` per contact into `Results_VTA_<strategy>E1C<i>/`.
# Superposition over current protocols + VTA metrics are computed
# afterwards by `evaluate_convergence_vta.py`.
#
# The refinement strategies mirror the plain VTA example
# (examples/ConvergenceStudy/VTA/.../run_convergence_study.py).  No
# neuron/pathway-derived netgen meshsize file is used here.
#
# Usage:
#   python run_convergence_study_vta.py                 # run all sequentially
#   python run_convergence_study_vta.py hp_refinement   # run one strategy
#   python run_convergence_study_vta.py VII VIII best    # mix names and numerals
#   python run_convergence_study_vta.py --list           # show available strategies

import argparse
import json
import logging
import os
from copy import deepcopy

import numpy as np

import ossdbs
from ossdbs.main import main_run

CONFIG_FILE = "oss-dbs_parameters_vta.json"


def remove_file_handler(logger):
    """Remove file handler from logger instance."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


def save_input_dict(base_input_dict):
    """Save the effective input dictionary alongside the results."""
    input_dict = deepcopy(base_input_dict)
    # add one layer so input files are found if re-run from the results dir
    input_dict["MaterialDistribution"]["MRIPath"] = os.path.join(
        "..", input_dict["MaterialDistribution"]["MRIPath"]
    )
    input_dict["StimSets"]["StimSetsFile"] = os.path.join(
        "..", input_dict["StimSets"]["StimSetsFile"]
    )
    input_dict["OutputPath"] = "./"
    with open(
        os.path.join(base_input_dict["OutputPath"], "input_dict.json"), "w"
    ) as fp:
        json.dump(input_dict, fp, indent=2)


def run_strategy(input_dict, logger):
    """Run a single convergence strategy (unit solves only, no PAM)."""
    # main_run dispatches to run_stim_sets because StimSets.Active is True.
    main_run(input_dict)
    save_input_dict(input_dict)
    remove_file_handler(logger)


def setup_base_config():
    """Load base configuration and precompute shared values."""
    with open(CONFIG_FILE) as fp:
        base = json.load(fp)

    electrode_name = base["Electrodes"][0]["Name"]

    # clean base state
    base["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    base["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
    base["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    base["Mesh"]["MaterialRefinementSteps"] = 0
    base["Mesh"]["HPRefinement"] = {"Active": False}
    _set_edge_size(base, 1e6)

    # load MRI for voxel sizes
    mri_image, _ = ossdbs.load_images(base)

    # edge sizes derived from electrode geometry
    lead_diameter = ossdbs.electrodes.default_electrode_parameters[
        electrode_name
    ].lead_diameter
    perimeter = np.pi * lead_diameter

    return base, mri_image, perimeter


def _set_edge_size(cfg, edge_size):
    """Set MaxMeshSizeEdge on all contacts.

    All contacts are activated in turn by the StimSets unit-solution
    loop, so edge refinement is applied to every contact (not just two).
    """
    for contact in cfg["Electrodes"][0]["Contacts"]:
        contact["MaxMeshSizeEdge"] = edge_size


# --- Strategy configuration functions ---
# Each returns a fully-configured, independent deepcopy of the base config.
# Definitions mirror the plain VTA convergence example.


def configure_default(base, mri_image, perimeter):
    """Default mesh."""
    cfg = deepcopy(base)
    cfg["OutputPath"] = "Results_VTA_default"
    return cfg


def configure_fine(base, mri_image, perimeter):
    """Fine mesh."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
    cfg["OutputPath"] = "Results_VTA_fine"
    return cfg


def configure_very_fine(base, mri_image, perimeter):
    """Very fine mesh."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
    cfg["OutputPath"] = "Results_VTA_very_fine"
    return cfg


def configure_material_refinement(base, mri_image, perimeter):
    """Very fine + 2x material refinement."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_VTA_material_refinement"
    return cfg


def configure_edge_refinement(base, mri_image, perimeter):
    """Default + edge refinement (perimeter/20)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 20.0)
    cfg["OutputPath"] = "Results_VTA_edge_refinement"
    return cfg


def configure_fine_edge_refinement(base, mri_image, perimeter):
    """Default + fine edge refinement (perimeter/50)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["OutputPath"] = "Results_VTA_fine_edge_refinement"
    return cfg


def configure_edge_voxel_refinement(base, mri_image, perimeter):
    """Fine edge (perimeter/50) + 10x voxel-size mesh limit."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 10.0 * min(mri_image.voxel_sizes)
    cfg["OutputPath"] = "Results_VTA_edge_voxel_refinement"
    return cfg


def configure_edge_single_material_refinement(base, mri_image, perimeter):
    """Fine edge (perimeter/50) + 1x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_VTA_edge_single_material_refinement"
    return cfg


def configure_edge_double_material_refinement(base, mri_image, perimeter):
    """Fine edge (perimeter/50) + 2x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_VTA_edge_double_material_refinement"
    return cfg


def configure_very_fine_edge_refinement(base, mri_image, perimeter):
    """Default + very fine edge refinement (perimeter/75)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["OutputPath"] = "Results_VTA_very_fine_edge_refinement"
    return cfg


def configure_fine_edge_single_material_refinement(base, mri_image, perimeter):
    """Very fine edge (perimeter/75) + 1x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_VTA_fine_edge_single_material_refinement"
    return cfg


def configure_fine_edge_double_material_refinement(base, mri_image, perimeter):
    """Very fine edge (perimeter/75) + 2x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_VTA_fine_edge_double_material_refinement"
    return cfg


def configure_hp_refinement(base, mri_image, perimeter):
    """Default mesh + HP refinement."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {"Active": True, "Levels": 2, "Factor": 0.125}
    cfg["OutputPath"] = "Results_VTA_hp_refinement"
    return cfg


def configure_hp_material_refinement(base, mri_image, perimeter):
    """Default mesh + HP ref. + 1x material ref."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {"Active": True, "Levels": 2, "Factor": 0.125}
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_VTA_hp_material_refinement"
    return cfg


def configure_best(base, mri_image, perimeter):
    """Gold standard: voxel mesh + material ref. + adaptive."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = min(mri_image.voxel_sizes)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    _set_edge_size(cfg, 1e6)
    cfg["Mesh"]["AdaptiveMeshRefinement"]["Active"] = True
    cfg["OutputPath"] = "Results_VTA_best"
    return cfg


# Strategy registry: (name, configure_function). Order mirrors the VTA example.
STRATEGIES = [
    ("default", configure_default),
    ("fine", configure_fine),
    ("very_fine", configure_very_fine),
    ("material_refinement", configure_material_refinement),
    ("edge_refinement", configure_edge_refinement),
    ("fine_edge_refinement", configure_fine_edge_refinement),
    ("edge_voxel_refinement", configure_edge_voxel_refinement),
    ("edge_single_material_refinement", configure_edge_single_material_refinement),
    ("edge_double_material_refinement", configure_edge_double_material_refinement),
    ("very_fine_edge_refinement", configure_very_fine_edge_refinement),
    (
        "fine_edge_single_material_refinement",
        configure_fine_edge_single_material_refinement,
    ),
    (
        "fine_edge_double_material_refinement",
        configure_fine_edge_double_material_refinement,
    ),
    ("hp_refinement", configure_hp_refinement),
    ("hp_material_refinement", configure_hp_material_refinement),
    ("best", configure_best),
]
STRATEGY_MAP = dict(STRATEGIES)

# Roman numeral lookup (I.. → strategy name)
_ROMAN = [
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
]
ROMAN_TO_NAME = {r: name for r, (name, _) in zip(_ROMAN, STRATEGIES, strict=False)}


def resolve_strategy_name(token):
    """Resolve a CLI token to a strategy name."""
    upper = token.upper()
    if upper in ROMAN_TO_NAME:
        return ROMAN_TO_NAME[upper]
    if token in STRATEGY_MAP:
        return token
    return None


def main():
    """Run selected convergence strategies sequentially."""
    parser = argparse.ArgumentParser(
        description="Run StimSets VTA convergence study strategies sequentially.",
    )
    parser.add_argument(
        "strategies",
        nargs="*",
        help="Strategies to run (by name or roman numeral). Default: all.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available strategies and exit.",
    )
    parser.add_argument(
        "--loglevel",
        type=int,
        default=logging.INFO,
        help="Logging verbosity (10=DEBUG, 20=INFO).",
    )
    args = parser.parse_args()

    if args.list:
        print("Available strategies:")
        for i, (name, fn) in enumerate(STRATEGIES):
            numeral = _ROMAN[i] if i < len(_ROMAN) else "?"
            print(f"  {numeral:>5s}  {name:45s}  {fn.__doc__}")
        return

    if args.strategies:
        to_run = []
        for token in args.strategies:
            resolved = resolve_strategy_name(token)
            if resolved is None:
                parser.error(f"Unknown strategy: {token!r}. Use --list to see options.")
            to_run.append(resolved)
    else:
        to_run = [name for name, _ in STRATEGIES]

    ossdbs.set_logger(level=args.loglevel)
    logger = logging.getLogger("ossdbs")
    base, mri_image, perimeter = setup_base_config()

    for name in to_run:
        print(f"\n{'=' * 60}")
        print(f"Running strategy: {name}")
        print(f"{'=' * 60}")
        cfg = STRATEGY_MAP[name](base, mri_image, perimeter)
        run_strategy(cfg, logger)

    print(f"\nCompleted {len(to_run)} strategy(ies).")


if __name__ == "__main__":
    main()
