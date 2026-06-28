# gold standard:
# impedance error 0.1% and maximum element size equal
# to mesh size

import json
import logging
import os
import sys

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


# Optional CLI: run only the named strategies (default: all). Names match the
# OutputPath suffix, e.g. "hp_double_material_refinement". Only the FEM solve
# is gated -- every base_input_dict mutation below still executes, so the
# cumulative mesh state is identical no matter which strategies are selected.
_SELECTED = {arg for arg in sys.argv[1:] if not arg.startswith("-")}


def run_selected(input_dict):
    """Run main_run only if this strategy was selected (or none were)."""
    name = input_dict["OutputPath"].replace("Results_VTA_", "")
    if _SELECTED and name not in _SELECTED:
        print(f"Skipping (not selected): {name}")
        return
    print(f"Running strategy: {name}")
    main_run(input_dict)
    remove_file_handler(_logger)


electrode_name = "Medtronic3389"

with open("../../base_settings.json") as fp:
    base_input_dict = json.load(fp)

# adjust pathes
base_input_dict["MaterialDistribution"]["MRIPath"] = os.path.join(
    "..", "..", base_input_dict["MaterialDistribution"]["MRIPath"]
)
# add electrode
base_input_dict["Electrodes"][0]["Name"] = electrode_name
# add contacts
contact_2_dict = BASE_CONTACT.copy()
contact_2_dict["Contact_ID"] = 2
contact_2_dict["Active"] = True
contact_2_dict["Voltage[V]"] = 1.0
base_input_dict["Electrodes"][0]["Contacts"].append(contact_2_dict)

# update lattice
base_input_dict["PointModel"]["Lattice"]["Shape"]["z"] = 90

# case grounding
surface_dict = {}
surface_dict["Name"] = "BrainSurface"
surface_dict["Active"] = True
surface_dict["Voltage[V]"] = 0.0
base_input_dict["Surfaces"] = [surface_dict]

# initially no adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"] = {}
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = False

# first refinement level: Default
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_default"
run_selected(base_input_dict)

# second refinement level: fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"
base_input_dict["OutputPath"] = "Results_VTA_fine"
run_selected(base_input_dict)

# third refinement level: very fine assumption
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "VeryFine"
base_input_dict["OutputPath"] = "Results_VTA_very_fine"
run_selected(base_input_dict)

# fourth refinement: material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["OutputPath"] = "Results_VTA_material_refinement"
run_selected(base_input_dict)
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
base_input_dict["OutputPath"] = "Results_VTA_edge_refinement"
run_selected(base_input_dict)

# sixth refinement: more edge refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
lead_diameter = ossdbs.electrodes.default_electrode_parameters[
    electrode_name
].lead_diameter
perimeter = np.pi * lead_diameter
edge_size = perimeter / 50.0
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_refinement"
run_selected(base_input_dict)

# eigth refinement: edge refinement + limit on voxel size
mri_image, _ = ossdbs.load_images(base_input_dict)
max_mesh_size = 10.0 * min(mri_image.voxel_sizes)
print(f"Imposing max mesh size of: {max_mesh_size:.2f}")
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
base_input_dict["OutputPath"] = "Results_VTA_edge_voxel_refinement"
run_selected(base_input_dict)
# undo mesh size imposing
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["OutputPath"] = "Results_VTA_edge_single_material_refinement"
run_selected(base_input_dict)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_edge_double_material_refinement"
run_selected(base_input_dict)
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
base_input_dict["OutputPath"] = "Results_VTA_very_fine_edge_refinement"
run_selected(base_input_dict)

# ninth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_single_material_refinement"
run_selected(base_input_dict)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# tenth refinement: material refinement + edge refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["OutputPath"] = "Results_VTA_fine_edge_double_material_refinement"
run_selected(base_input_dict)
# reset material refinement
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0

# HP refinement: default mesh + HP refinement
base_input_dict["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = 1e6
base_input_dict["Mesh"]["HPRefinement"] = {
    "Active": True,
    "Levels": 2,
    "Factor": 0.125,
}
base_input_dict["OutputPath"] = "Results_VTA_hp_refinement"
run_selected(base_input_dict)

# HP + material refinement: default mesh + HP ref. + 1x material ref.
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
base_input_dict["OutputPath"] = "Results_VTA_hp_material_refinement"
run_selected(base_input_dict)

# HP + double material refinement: default mesh + HP ref. + 2x material ref.
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 2
base_input_dict["OutputPath"] = "Results_VTA_hp_double_material_refinement"
run_selected(base_input_dict)
# reset
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 0
base_input_dict["Mesh"]["HPRefinement"] = {"Active": False}

# finest level: voxel size + adaptive mesh refinement
max_mesh_size = min(mri_image.voxel_sizes)
base_input_dict["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max_mesh_size
# one material refinement step
base_input_dict["Mesh"]["MaterialRefinementSteps"] = 1
# reset edge size
edge_size = 1e6
base_input_dict["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size
# adaptive refinement
base_input_dict["Mesh"]["AdaptiveMeshRefinement"]["Active"] = True

base_input_dict["OutputPath"] = "Results_VTA_best"
run_selected(base_input_dict)
