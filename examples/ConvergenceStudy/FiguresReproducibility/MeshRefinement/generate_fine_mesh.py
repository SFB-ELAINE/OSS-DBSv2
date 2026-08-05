# Generate the global "fine" / "very_fine" meshes for the mesh-refinement
# figure (a visible, uniform refinement to contrast with the default mesh).

import json
import logging

import ngsolve

import ossdbs

ossdbs.set_logger(logging.INFO)

electrode_name = "Medtronic3389"

with open("oss_dbs_parameters.json") as fp:
    base_cfg = json.load(fp)
base_cfg["Electrodes"][0]["Name"] = electrode_name
base_cfg["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}

mri_image, _ = ossdbs.load_images(base_cfg)
electrodes = ossdbs.api.generate_electrodes(base_cfg)
brain_model = ossdbs.api.build_brain_model(base_cfg, mri_image)
try:
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)
except RuntimeError:
    brain_model = ossdbs.api.build_brain_model(
        base_cfg, mri_image, rotate_initial_geo=True
    )
    geometry = ossdbs.ModelGeometry(brain_model, electrodes)

solver = ossdbs.api.prepare_solver(base_cfg)
dielectric_properties = ossdbs.api.prepare_dielectric_properties(base_cfg)
materials = base_cfg["MaterialDistribution"]["MRIMapping"]
conductivity = ossdbs.ConductivityCF(
    mri_image,
    brain_model.brain_region,
    dielectric_properties,
    materials,
    geometry.encapsulation_layers,
    complex_data=base_cfg["EQSMode"],
)

for hypothesis, out in [
    ("Fine", "Results_PAM_fine"),
    ("VeryFine", "Results_PAM_very_fine"),
]:
    cfg = json.loads(json.dumps(base_cfg))
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = hypothesis
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    cfg["OutputPath"] = out
    with ngsolve.TaskManager():
        volume_conductor = ossdbs.api.prepare_volume_conductor_model(
            cfg, geometry, conductivity, solver
        )
        print(f"{out}: {volume_conductor.mesh.ngsolvemesh.ne} elements")
        volume_conductor.export_material_distribution_to_vtk()
