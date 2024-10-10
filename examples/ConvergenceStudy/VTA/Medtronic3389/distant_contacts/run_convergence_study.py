# gold standard:
# impedance error 0.1% and maximum element size equal
# to mesh size

import json
import logging

import numpy as np

import ossdbs
from ossdbs.main import main_run


def remove_file_handler(logger):
    """Remove file handler from logger instance."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


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

with open("../../base_settings.json") as fp:
    base_input_dict = json.load(fp)

# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name
# add contacts
contact_1_dict = BASE_CONTACT.copy()
contact_1_dict["Contact_ID"] = 1
contact_1_dict["Active"] = True
contact_1_dict["Voltage[V]"] = 1.0
contact_2_dict = BASE_CONTACT.copy()
contact_2_dict["Contact_ID"] = 4
contact_2_dict["Active"] = True
base_input_dict["Electrodes"][0]["Contacts"].append(contact_1_dict)
base_input_dict["Electrodes"][0]["Contacts"].append(contact_2_dict)

# update lattice
base_input_dict["PointModel"]["Lattice"]["Shape"]["z"] = 90

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# first refinement level: Default
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_default"
main_run(base_input_dict)
remove_file_handler(_logger)

# second refinement level: fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
base_input_dict["OutputPath"] = "Results_VTA_fine"
main_run(base_input_dict)
remove_file_handler(_logger)

# third refinement level: very fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
base_input_dict["OutputPath"] = "Results_VTA_very_fine"
main_run(base_input_dict)
remove_file_handler(_logger)

# fourth refinement: material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["OutputPath"] = "Results_VTA_material_refinement"
main_run(base_input_dict)
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
base_input_dict["OutputPath"] = "Results_VTA_edge_refinement"
main_run(base_input_dict)
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
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)

# eigth refinement: edge refinement + limit on voxel size
mri_image, _ = ossdbs.load_images(base_input_dict)
max_mesh_size = 10.0 * min(mri_image.voxel_sizes)
print(f"Imposing max mesh size of: {max_mesh_size:.2f}")
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
base_input_dict["OutputPath"] = "Results_VTA_edge_voxel_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_edge_single_material_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_edge_double_material_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# seventh refinement: very fine edge refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 75.0
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
base_input_dict["OutputPath"] = "Results_VTA_very_fine_edge_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_single_material_refinement"
main_run(base_input_dict)
remove_file_handler(_logger)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# tenth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_double_material_refinement"
main_run(base_input_dict)
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

base_input_dict["OutputPath"] = "Results_VTA_best"
main_run(base_input_dict)
