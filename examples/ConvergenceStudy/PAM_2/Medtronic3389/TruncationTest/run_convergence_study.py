# gold standard:
# impedance error 0.1% and maximum element size equal
# to mesh size

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
        os.path.join(base_input_dict["OutputPath"], "input_dict.json"), "w"
    ) as fp:
        json.dump(input_dict, fp, indent=2)


ossdbs.set_logger()
_logger = logging.getLogger("ossdbs")

electrode_name = "BostonScientificVerciseDirected"

with open("../../oss-dbs_parameters.json") as fp:
    base_input_dict = json.load(fp)


# adjust pathes
base_input_dict["MaterialDistribution"]["MRIPath"] = os.path.join(
    "..", "..", base_input_dict["MaterialDistribution"]["MRIPath"]
)
base_input_dict["PointModel"]["Pathway"]["FileName"] = os.path.join(
    "..", "..", base_input_dict["PointModel"]["Pathway"]["FileName"]
)
base_input_dict["PathwayFile"] = os.path.join(
    "..", "..", base_input_dict["PathwayFile"]
)

# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name

# for PAM
base_input_dict["Scaling"] = 1.0
base_input_dict["ScalingIndex"] = None

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# refinement of edges and material
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 50.0
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size

base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_PAM_no_truncation"
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

for truncation_ratio in [5, 10, 20, 30]:
    base_input_dict["TruncateAfterActivePartRatio"] = float(truncation_ratio)
    base_input_dict["OutputPath"] = f"Results_PAM_truncation_{truncation_ratio}"
    main_run(base_input_dict)
    save_input_dict(base_input_dict)
    ossdbs.api.run_PAM(base_input_dict)
    remove_file_handler(_logger)
