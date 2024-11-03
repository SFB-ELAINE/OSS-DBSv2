# gold standard:
# impedance error 0.1% and maximum element size equal
# to mesh size

import json
import logging
import os
from copy import deepcopy

import ngsolve
import numpy as np

import ossdbs


def remove_file_handler(logger):
    """Remove file handler from logger instance."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


def save_input_dict(base_input_dict):
    """Save input dictionary for convergence run."""
    input_dict = deepcopy(base_input_dict)
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

ossdbs.set_logger(logging.INFO)
_logger = logging.getLogger("ossdbs")

electrode_name = "Medtronic3389"

with open("oss_dbs_parameters.json") as fp:
    base_input_dict = json.load(fp)


# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name

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
electrodes = ossdbs.api.generate_electrodes(base_input_dict)
brain_model = ossdbs.api.build_brain_model(base_input_dict, mri_image)
try:
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)
except RuntimeError:
    brain_model = ossdbs.api.build_brain_model(
        base_input_dict, mri_image, rotate_initial_geo=True
    )
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)


solver = ossdbs.api.prepare_solver(base_input_dict)
dielectric_properties = ossdbs.api.prepare_dielectric_properties(base_input_dict)
materials = base_input_dict["MaterialDistribution"]["MRIMapping"]
conductivity = ossdbs.ConductivityCF(
    mri_image,
    brain_model.brain_region,
    dielectric_properties,
    materials,
    geometry.encapsulation_layers,
    complex_data=base_input_dict["EQSMode"],
)


with ngsolve.TaskManager():
    volume_conductor = ossdbs.api.prepare_volume_conductor_model(
        base_input_dict, geometry, conductivity, solver
    )
    print(volume_conductor.mesh.ngsolvemesh.ne)
    print("Save material")
    volume_conductor.export_material_distribution_to_vtk()

# refinement level: default and refinement around pathway
base_input_dict["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = "meshsizes.txt"
base_input_dict["OutputPath"] = "Results_PAM_default_meshsize"
with ngsolve.TaskManager():
    volume_conductor = ossdbs.api.prepare_volume_conductor_model(
        base_input_dict, geometry, conductivity, solver
    )
    print(volume_conductor.mesh.ngsolvemesh.ne)
    print("Save material")
    volume_conductor.export_material_distribution_to_vtk()
# reset mesh size file name
base_input_dict["Mesh"]["MeshingHypothesis"]["MeshSizeFilename"] = ""


# fourth refinement: material refinement
base_input_dict["OutputPath"] = "Results_PAM_material_refinement"
with ngsolve.TaskManager():
    volume_conductor = ossdbs.api.prepare_volume_conductor_model(
        base_input_dict, geometry, conductivity, solver
    )
    # two refinement steps
    volume_conductor.refine_mesh_by_material(2)
    print(volume_conductor.mesh.ngsolvemesh.ne)
    print("Save material")
    volume_conductor.export_material_distribution_to_vtk()

# sixth refinement: more edge refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 50.0
print(edge_size)
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["Electrodes"][0]["Contacts"][1]["MaxMeshSizeEdge"] = edge_size
base_input_dict["OutputPath"] = "Results_PAM_fine_edge_refinement"
ossdbs.api.set_contact_and_encapsulation_layer_properties(base_input_dict, geometry)


with ngsolve.TaskManager():
    volume_conductor = ossdbs.api.prepare_volume_conductor_model(
        base_input_dict, geometry, conductivity, solver
    )
    geo = volume_conductor.mesh.ngsolvemesh.ngmesh.GetGeometry()
    for edge in geo.shape.edges:
        print(edge.name, edge.maxh)
    print(volume_conductor.mesh.ngsolvemesh.ne)
    print("Save material")
    volume_conductor.export_material_distribution_to_vtk()
