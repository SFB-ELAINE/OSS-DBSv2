# Convergence study: signal truncation strategies for PAM
#
# Usage:
#   python run_convergence_study.py                    # run all sequentially

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
    """Load base configuration and precompute shared values."""
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
# No strategy depends on state left by a previous one.
def configure_hp_material_refinement(base, mri_image, perimeter):
    """Default mesh + HP ref. + 1x material ref."""
    cfg = deepcopy(base)
    cfg["Mesh"]["HPRefinement"] = {
        "Active": True,
        "Levels": 2,
        "Factor": 0.125,
    }
    cfg["Mesh"]["MaterialRefinementSteps"] = 1
    cfg["OutputPath"] = "Results_PAM_hp_material_refinement"
    return cfg


def main():
    """Run sequentially with different truncation ratios."""
    parser = argparse.ArgumentParser(
        description=("Run PAM convergence study on truncation ratio."),
    )
    parser.add_argument(
        "--loglevel",
        type=int,
        default=logging.INFO,
        help="Logging verbosity (10=DEBUG, 20=INFO).",
    )
    args = parser.parse_args()

    # Truncation ratios to be checked (None = no truncation)
    to_run = [None, 5, 10, 20, 30]

    ossdbs.set_logger(level=args.loglevel)
    logger = logging.getLogger("ossdbs")
    base, mri_image, perimeter = setup_base_config()

    for ratio in to_run:
        label = "none" if ratio is None else str(ratio)
        print(f"\n{'=' * 60}")
        print(f"Running strategy: truncation_{label}")
        print(f"{'=' * 60}")
        cfg = configure_hp_material_refinement(base, mri_image, perimeter)
        if ratio is not None:
            cfg["TruncateAfterActivePartRatio"] = float(ratio)
        cfg["OutputPath"] = f"Results_PAM_truncation_{label}"
        run_strategy(cfg, logger)

    print(f"\nCompleted {len(to_run)} strategy(ies).")


if __name__ == "__main__":
    main()
