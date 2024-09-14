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
        "..", input_dict["MaterialDistribution"]["MRIPath"]
    )
    input_dict["PointModel"]["Pathway"]["FileName"] = os.path.join(
        "..", input_dict["PointModel"]["Pathway"]["FileName"]
    )
    input_dict["PathwayFile"] = os.path.join("..", input_dict["PathwayFile"])
    input_dict["OutputPath"] = "./"

    with open(
        os.path.join(base_input_dict["OutputPath"], "input_dict.json"), "w"
    ) as fp:
        json.dump(input_dict, fp, indent=2)


BASE_CONTACT = {
    "Contact_ID": 1,
    "Active": False,
    "Current[A]": 0.0,
    "Voltage[V]": 0.0,
    "Floating": False,
}

ossdbs.set_logger()
_logger = logging.getLogger("ossdbs")

electrode_name = "Medtronic3389"

with open("../oss_dbs_parameters.json") as fp:
    base_input_dict = json.load(fp)


# adjust pathes
base_input_dict["MaterialDistribution"]["MRIPath"] = os.path.join(
    "..", base_input_dict["MaterialDistribution"]["MRIPath"]
)
base_input_dict["PointModel"]["Pathway"]["FileName"] = os.path.join(
    "..", base_input_dict["PointModel"]["Pathway"]["FileName"]
)
base_input_dict["PathwayFile"] = os.path.join("..", base_input_dict["PathwayFile"])

# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name

# for PAM
base_input_dict["Scaling"] = 1.0
base_input_dict["ScalingIndex"] = None

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# first refinement level: Default
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_PAM_default"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# second refinement level: fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
base_input_dict["OutputPath"] = "Results_PAM_fine"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# third refinement level: very fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
base_input_dict["OutputPath"] = "Results_PAM_very_fine"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# fourth refinement: material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["OutputPath"] = "Results_PAM_material_refinement"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# fifth refinement: edge refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 20.0
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
base_input_dict["OutputPath"] = "Results_PAM_edge_refinement"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# sixth refinement: more edge refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 50.0
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
base_input_dict["OutputPath"] = "Results_PAM_fine_edge_refinement"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# seventh refinement: edge refinement + limit on voxel size
mri_image, _ = ossdbs.load_images(base_input_dict)
max_mesh_size = 10.0 * min(mri_image.voxel_sizes)
print(f"Imposing max mesh size of: {max_mesh_size:.2f}")
base_input_dict["OutputPath"] = "Results_PAM_edge_voxel_refinement"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# eigth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_PAM_edge_material_refinement"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# finest level: voxel size + adaptive mesh refinement
max_mesh_size = min(mri_image.voxel_sizes)
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
# one material refinement step
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
# reset edge size
edge_size = 1e6
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
# adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = True
base_input_dict["OutputPath"] = "Results_PAM_best"
save_input_dict(base_input_dict)
main_run(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
