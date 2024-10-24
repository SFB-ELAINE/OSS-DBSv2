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


ossdbs.set_logger()
_logger = logging.getLogger("ossdbs")

with open("oss-dbs_parameters.json") as fp:
    base_input_dict = json.load(fp)

electrode_name = base_input_dict["Electrodes"][0]["Name"]

# do not scale solution
base_input_dict["Scaling"] = 1.0

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# mri image
mri_image, _ = ossdbs.load_images(base_input_dict)
# pathway info
pw = ossdbs.point_analysis.Pathway(base_input_dict["PointModel"]["Pathway"]["FileName"])
pw.write_netgen_meshsize_file(
    meshsize=min(mri_image.voxel_sizes), filename="meshsizes.txt"
)

# first refinement level: Default
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_PAM_default"
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# refinement level: default and refinement around pathway
base_input_dict["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = "meshsizes.txt"
base_input_dict["OutputPath"] = "Results_PAM_default_meshsize"
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)
# reset mesh size file name
base_input_dict["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = ""

# second refinement level: fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
base_input_dict["OutputPath"] = "Results_PAM_fine"
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# fourth refinement: material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["OutputPath"] = "Results_PAM_material_refinement"
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

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
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
remove_file_handler(_logger)

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_PAM_edge_single_material_refinement"
main_run(base_input_dict)
save_input_dict(base_input_dict)
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
main_run(base_input_dict)
save_input_dict(base_input_dict)
ossdbs.api.run_PAM(base_input_dict)
