# Convergence study: nonfloating (voltage-controlled) reference
#
# Same mesh strategies as run_convergence_study.py but with
# voltage-controlled contacts (Dirichlet BCs) instead of
# current-controlled floating contacts. This uses the nonfloating
# volume conductor path which is unaffected by the 1e-3 area
# scaling issue in floating.py.
#
# Contact setup:
#   E1C2/3/4 =  1 V (cathodes, active)
#   E1C8     =  0 V (anode/ground, active)
#   Others   = floating, no current
#
# Usage:
#   python run_convergence_study_nonfloating.py
#   python run_convergence_study_nonfloating.py hp_refinement
#   python run_convergence_study_nonfloating.py --list

import argparse
import json
import logging
import os
from copy import deepcopy

import numpy as np

import ossdbs
from ossdbs.main import main_run


def remove_file_handler(logger):
    """Remove file handler from logger instance."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


def save_input_dict(base_input_dict):
    """Save input dictionary for convergence run."""
    input_dict = deepcopy(base_input_dict)
    # add one layer to make sure input files get found
    input_dict["MaterialDistribution"]["MRIPath"] = os.path.join(
        "..", "..", input_dict["MaterialDistribution"]["MRIPath"]
    )
    input_dict["PointModel"]["Pathway"]["FileName"] = os.path.join(
        "..", "..", input_dict["PointModel"]["Pathway"]["FileName"]
    )
    input_dict["PathwayFile"] = os.path.join("..", "..", input_dict["PathwayFile"])
    input_dict["OutputPath"] = "./"

    with open(
        os.path.join(base_input_dict["OutputPath"], "input_dict.json"),
        "w",
    ) as fp:
        json.dump(input_dict, fp, indent=2)


def run_strategy(input_dict, logger):
    """Run a single convergence strategy."""
    main_run(input_dict)
    save_input_dict(input_dict)
    ossdbs.api.run_PAM(input_dict)
    remove_file_handler(logger)


def setup_base_config():
    """Load base configuration and convert to voltage-controlled."""
    electrode_name = "BostonScientificVerciseDirected"

    with open("../../oss-dbs_parameters.json") as fp:
        base = json.load(fp)

    # adjust paths
    base["MaterialDistribution"]["MRIPath"] = os.path.join(
        "..", "..", base["MaterialDistribution"]["MRIPath"]
    )
    base["PointModel"]["Pathway"]["FileName"] = os.path.join(
        "..", "..", base["PointModel"]["Pathway"]["FileName"]
    )
    base["PathwayFile"] = os.path.join("..", "..", base["PathwayFile"])

    # for PAM
    base["Scaling"] = 1.0
    base["ScalingIndex"] = None

    # --- Convert to voltage-controlled (nonfloating) ---
    base["StimulationSignal"]["CurrentControlled"] = False
    for contact in base["Electrodes"][0]["Contacts"]:
        cid = contact["Contact_ID"]
        contact["Current[A]"] = 0.0
        if cid in (2, 3, 4):
            # cathodes: -1 V
            contact["Active"] = True
            contact["Floating"] = False
            contact["Voltage[V]"] = 3.5
        elif cid == 8:
            # anode / ground: 0 V
            contact["Active"] = True
            contact["Floating"] = False
            contact["Voltage[V]"] = 0.0
        else:
            # inactive contacts float freely
            contact["Active"] = False
            contact["Floating"] = True
            contact["Voltage[V]"] = 0.0

    # clean base state
    base["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    base["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
    base["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    base["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = ""
    base["Mesh"]["MaterialRefinementSteps"] = 0

    # load MRI for voxel sizes
    mri_image, _ = ossdbs.load_images(base)

    # pathway meshsize file (reuse if already generated)
    if not os.path.isfile("meshsizes.txt"):
        pw = ossdbs.point_analysis.Pathway(base["PointModel"]["Pathway"]["FileName"])
        pw.write_netgen_meshsize_file(
            meshsize=min(mri_image.voxel_sizes),
            filename="meshsizes.txt",
        )

    # compute edge sizes from electrode geometry
    lead_diameter = ossdbs.electrodes.default_electrode_parameters[
        electrode_name
    ].lead_diameter
    perimeter = np.pi * lead_diameter

    return base, mri_image, perimeter


def _set_edge_size(cfg, edge_size):
    """Set MaxMeshSizeEdge on the two active contacts."""
    cfg["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
    cfg["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size


# --- Strategy configuration functions ---
# Each returns a fully-configured deepcopy of the base config.
# Output paths use NF_ prefix to distinguish from floating results.


def configure_default(base, mri_image, perimeter):
    """Default mesh."""
    cfg = deepcopy(base)
    cfg["OutputPath"] = "Results_NF_default"
    return cfg


def configure_default_meshsize(base, mri_image, perimeter):
    """Default + pathway meshsize."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = "meshsizes.txt"
    cfg["OutputPath"] = "Results_NF_default_meshsize"
    return cfg


def configure_fine(base, mri_image, perimeter):
    """Fine mesh."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
    cfg["OutputPath"] = "Results_NF_fine"
    return cfg


def configure_very_fine(base, mri_image, perimeter):
    """VeryFine mesh."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
    cfg["OutputPath"] = "Results_NF_very_fine"
    return cfg


def configure_material_refinement(base, mri_image, perimeter):
    """VeryFine + 2x material refinement."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_NF_material_refinement"
    return cfg


def configure_edge_refinement(base, mri_image, perimeter):
    """Edge refinement (perimeter/20)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 20.0)
    cfg["OutputPath"] = "Results_NF_edge_refinement"
    return cfg


def configure_fine_edge_refinement(base, mri_image, perimeter):
    """Fine edge refinement (perimeter/50)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["OutputPath"] = "Results_NF_fine_edge_refinement"
    return cfg


def configure_edge_voxel_refinement(base, mri_image, perimeter):
    """Fine edge + voxel-based max mesh size."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    max_mesh_size = 10.0 * min(mri_image.voxel_sizes)
    print(f"Imposing max mesh size of: {max_mesh_size:.2f}")
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
    cfg["OutputPath"] = "Results_NF_edge_voxel_refinement"
    return cfg


def configure_edge_single_material_refinement(base, mri_image, perimeter):
    """Fine edge + 1x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_NF_edge_single_material_refinement"
    return cfg


def configure_edge_double_material_refinement(base, mri_image, perimeter):
    """Fine edge + 2x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 50.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_NF_edge_double_material_refinement"
    return cfg


def configure_very_fine_edge_refinement(base, mri_image, perimeter):
    """Very fine edge refinement (perimeter/75)."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["OutputPath"] = "Results_NF_very_fine_edge_refinement"
    return cfg


def configure_edge_meshsize(base, mri_image, perimeter):
    """Very fine edge + pathway meshsize."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = "meshsizes.txt"
    cfg["OutputPath"] = "Results_NF_edge_meshsize"
    return cfg


def configure_fine_edge_single_material_refinement(base, mri_image, perimeter):
    """Very fine edge + 1x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_NF_fine_edge_single_material_refinement"
    return cfg


def configure_fine_edge_double_material_refinement(base, mri_image, perimeter):
    """Very fine edge + 2x material refinement."""
    cfg = deepcopy(base)
    _set_edge_size(cfg, perimeter / 75.0)
    cfg["Mesh"]["MaterialRefinementSteps"] = 2
    cfg["OutputPath"] = "Results_NF_fine_edge_double_material_refinement"
    return cfg


def configure_hp_refinement(base, mri_image, perimeter):
    """Default mesh + HP refinement."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {
        "Active": True,
        "Levels": 2,
        "Factor": 0.125,
    }
    cfg["OutputPath"] = "Results_NF_hp_refinement"
    return cfg


def configure_hp_material_refinement(base, mri_image, perimeter):
    """Default mesh + HP ref. + 1x material ref."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {
        "Active": True,
        "Levels": 2,
        "Factor": 0.125,
    }
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_NF_hp_material_refinement"
    return cfg


def configure_best(base, mri_image, perimeter):
    """Gold standard: voxel mesh + material ref. + adaptive."""
    cfg = deepcopy(base)
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = min(mri_image.voxel_sizes)
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    _set_edge_size(cfg, 1e6)
    cfg["Mesh"]["AdaptiveMeshRefinement"]["Active"] = True
    cfg["OutputPath"] = "Results_NF_best"
    return cfg


# Strategy registry
STRATEGIES = [
    ("default", configure_default),
    ("default_meshsize", configure_default_meshsize),
    ("fine", configure_fine),
    ("very_fine", configure_very_fine),
    ("material_refinement", configure_material_refinement),
    ("edge_refinement", configure_edge_refinement),
    ("fine_edge_refinement", configure_fine_edge_refinement),
    ("edge_voxel_refinement", configure_edge_voxel_refinement),
    (
        "edge_single_material_refinement",
        configure_edge_single_material_refinement,
    ),
    (
        "edge_double_material_refinement",
        configure_edge_double_material_refinement,
    ),
    (
        "very_fine_edge_refinement",
        configure_very_fine_edge_refinement,
    ),
    ("edge_meshsize", configure_edge_meshsize),
    (
        "fine_edge_single_material_refinement",
        configure_fine_edge_single_material_refinement,
    ),
    (
        "fine_edge_double_material_refinement",
        configure_fine_edge_double_material_refinement,
    ),
    ("hp_refinement", configure_hp_refinement),
    (
        "hp_material_refinement",
        configure_hp_material_refinement,
    ),
    ("best", configure_best),
]
STRATEGY_MAP = dict(STRATEGIES)

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
    "XVI",
    "XVII",
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
        description=(
            "Run nonfloating (voltage-controlled) PAM "
            "convergence study strategies sequentially."
        ),
    )
    parser.add_argument(
        "strategies",
        nargs="*",
        help=("Strategies to run (by name or roman numeral). Default: all."),
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

    # Resolve and validate requested strategies
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
