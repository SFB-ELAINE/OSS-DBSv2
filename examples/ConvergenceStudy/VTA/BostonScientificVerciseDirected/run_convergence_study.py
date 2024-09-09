# gold standard:
# impedance error 0.1% and maximum element size equal
# to mesh size

import json

import numpy as np

import ossdbs
from ossdbs.main import main_run

BASE_CONTACT = {
    "Contact_ID": 1,
    "Active": False,
    "Current[A]": 0.0,
    "Voltage[V]": 0.0,
    "Floating": False,
}

ossdbs.set_logger()
electrode_name = "BostonScientificVerciseDirected"

with open("../base_settings.json") as fp:
    base_input_dict = json.load(fp)

# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name
# add contacts
contact_1_dict = BASE_CONTACT.copy()
contact_1_dict["Contact_ID"] = 1
contact_1_dict["Active"] = True
contact_1_dict["Voltage[V]"] = 1.0
contact_2_dict = BASE_CONTACT.copy()
contact_2_dict["Contact_ID"] = 2
contact_2_dict["Active"] = True
base_input_dict["Electrodes"][0]["Contacts"].append(contact_1_dict)
base_input_dict["Electrodes"][0]["Contacts"].append(contact_2_dict)

# update grid
base_input_dict["PointModel"]["Lattice"]["Shape"]["x"] = 60
base_input_dict["PointModel"]["Lattice"]["Shape"]["y"] = 60
base_input_dict["PointModel"]["Lattice"]["Shape"]["z"] = 120
base_input_dict["PointModel"]["Lattice"]["PointDistance[mm]"] = 0.125

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# first refinement level: Default
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_default"
main_run(base_input_dict)

# second refinement level: fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
base_input_dict["OutputPath"] = "Results_VTA_fine"
main_run(base_input_dict)

# third refinement level: very fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
base_input_dict["OutputPath"] = "Results_VTA_very_fine"
main_run(base_input_dict)

# fourth refinement: edge refinement
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

# fifth refinement: more edge refinement
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

# sixth refinement: edge refinement + limit on voxel size
mri_image, _ = ossdbs.load_images(base_input_dict)
max_mesh_size = 10.0 * min(mri_image.voxel_sizes)
print(f"Imposing max mesh size of: {max_mesh_size:.2f}")
base_input_dict["OutputPath"] = "Results_VTA_edge_voxel_refinement"
main_run(base_input_dict)

# finest level: voxel size + adaptive mesh refinement
max_mesh_size = min(mri_image.voxel_sizes)
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
# reset edge size
edge_size = 1e6
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
# adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = True

base_input_dict["OutputPath"] = "Results_VTA_best"
main_run(base_input_dict)
