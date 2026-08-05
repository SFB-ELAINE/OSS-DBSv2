# Generate the hp-refined mesh for the mesh-refinement figure.
# HP refinement must be active *before* the electrode is built so the
# contact edges are marked as singular for RefineHP, so this does its own
# geometry setup (generate_meshes.py reuses one geometry for all strategies).

import json
import logging

import ngsolve

import ossdbs

ossdbs.set_logger(logging.INFO)

electrode_name = "Medtronic3389"

with open("oss_dbs_parameters.json") as fp:
    cfg = json.load(fp)

cfg["Electrodes"][0]["Name"] = electrode_name
cfg["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
cfg["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
# HP refinement active from the start so the contacts are marked singular
cfg["Mesh"]["HPRefinement"] = {"Active": True, "Levels": 2, "Factor": 0.125}
cfg["OutputPath"] = "Results_PAM_hp_refinement"

mri_image, _ = ossdbs.load_images(cfg)
electrodes = ossdbs.api.generate_electrodes(cfg)
brain_model = ossdbs.api.build_brain_model(cfg, mri_image)
try:
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)
except RuntimeError:
    brain_model = ossdbs.api.build_brain_model(cfg, mri_image, rotate_initial_geo=True)
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)

solver = ossdbs.api.prepare_solver(cfg)
dielectric_properties = ossdbs.api.prepare_dielectric_properties(cfg)
materials = cfg["MaterialDistribution"]["MRIMapping"]
conductivity = ossdbs.ConductivityCF(
    mri_image,
    brain_model.brain_region,
    dielectric_properties,
    materials,
    geometry.encapsulation_layers,
    complex_data=cfg["EQSMode"],
)

with ngsolve.TaskManager():
    volume_conductor = ossdbs.api.prepare_volume_conductor_model(
        cfg, geometry, conductivity, solver
    )
    # deferred HP refinement (after any bisection-based refinement)
    volume_conductor.mesh.apply_hp_refinement()
    print("elements:", volume_conductor.mesh.ngsolvemesh.ne)
    print("Save material")
    volume_conductor.export_material_distribution_to_vtk()
