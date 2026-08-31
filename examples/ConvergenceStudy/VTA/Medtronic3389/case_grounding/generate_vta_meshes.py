"""Generate the VTA meshes for the mesh-refinement grid figure.

Mirrors ``FiguresReproducibility/MeshRefinement/generate_meshes.py`` but for
the VTA/Medtronic3389 case-grounding scenario. For each strategy it builds the
volume-conductor mesh and writes the *exact* ngsolve mesh + material labels to
``Results_VTA_<strategy>/material.vtu`` via
``export_material_distribution_to_vtk`` (subdivision 0 -> one VTK cell per
element, so the wireframe is the real FEM mesh).

It does NOT solve the FEM problem and does NOT run the VTA point analysis or
``ExportVTK``; only the material-distribution mesh needed by
``render_vta_mesh_grid_svg.py`` is produced.

The "best" panel here is the finest *static* mesh (voxel-size limit + material
refinement); the study's true ``Results_VTA_best`` adds a solution-driven
adaptive step, which would require a full solve.

Run from this directory (Medtronic3389/case_grounding); strategies run strictly
sequentially, one after the other.
"""

import gc
import json
import logging

import ngsolve
import numpy as np

import ossdbs
from ossdbs.utils.settings import Settings

ossdbs.set_logger(logging.INFO)
_logger = logging.getLogger("ossdbs")

ELECTRODE_NAME = "Medtronic3389"

BASE_CONTACT = {
    "Contact_ID": 1,
    "Active": False,
    "Current[A]": 0.0,
    "Voltage[V]": 0.0,
    "Floating": False,
}


def base_config():
    """Load base_settings and replicate the case-grounding study setup."""
    with open("../../base_settings.json") as fp:
        cfg = json.load(fp)
    cfg["MaterialDistribution"]["MRIPath"] = "../../segmask.nii.gz"
    cfg["Electrodes"][0]["Name"] = ELECTRODE_NAME
    # active contact (contact 2), as in run_convergence_study.py
    contact = BASE_CONTACT.copy()
    contact["Contact_ID"] = 2
    contact["Active"] = True
    contact["Voltage[V]"] = 1.0
    cfg["Electrodes"][0]["Contacts"].append(contact)
    # case grounding
    cfg["Surfaces"] = [{"Name": "BrainSurface", "Active": True, "Voltage[V]": 0.0}]
    # no adaptive / solution-driven refinement for pure meshing
    cfg["Mesh"]["AdaptiveMeshRefinement"] = {"Active": False}
    cfg["Mesh"]["MeshingHypothesis"]["Type"] = "Default"
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = 1e6
    # fill in defaults (DielectricModel.CustomParameters, etc.) like main_run
    return Settings(cfg).complete_settings()


def build_geometry(cfg, mri_image):
    """Fresh electrodes + brain + model geometry (with rotate fallback).

    A fresh geometry per strategy keeps contact/hp mesh-size properties from
    leaking between strategies (mesh sizes are snapshotted into OCCGeometry).
    """
    electrodes = ossdbs.api.generate_electrodes(cfg)
    brain_model = ossdbs.api.build_brain_model(cfg, mri_image)
    try:
        geometry = ossdbs.ModelGeometry(brain_model, electrodes)
    except RuntimeError:
        brain_model = ossdbs.api.build_brain_model(
            cfg, mri_image, rotate_initial_geo=True
        )
        geometry = ossdbs.ModelGeometry(brain_model, electrodes)
    # apply contact/edge mesh sizes onto the geometry before OCCGeometry is
    # built (ModelGeometry builds it lazily on first .geometry access)
    ossdbs.api.set_contact_and_encapsulation_layer_properties(cfg, geometry)
    return geometry, brain_model


def export_strategy(cfg, mri_image, material_steps=0, hp=False):
    """Build one strategy's mesh (no solve) and write its material.vtu."""
    geometry, brain_model = build_geometry(cfg, mri_image)
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
        vc = ossdbs.api.prepare_volume_conductor_model(
            cfg, geometry, conductivity, solver
        )
        if material_steps:
            vc.refine_mesh_by_material(material_steps)
        if hp:
            vc.mesh.apply_hp_refinement()
        _logger.info(f"{cfg['OutputPath']}: {vc.mesh.ngsolvemesh.ne} elements")
        vc.export_material_distribution_to_vtk()
    del vc, conductivity, geometry, brain_model, solver
    gc.collect()


def main():
    """Generate all six panels' meshes sequentially."""
    lead_diameter = ossdbs.electrodes.default_electrode_parameters[
        ELECTRODE_NAME
    ].lead_diameter
    edge_size = np.pi * lead_diameter / 20.0

    cfg = base_config()
    mri_image, _ = ossdbs.load_images(cfg)
    voxel_size = float(min(mri_image.voxel_sizes))

    # (label, output dir, mutator) -- run strictly sequentially
    def default(c):
        pass

    def fine(c):
        c["Mesh"]["MeshingHypothesis"]["Type"] = "Fine"

    def edge(c):
        c["Electrodes"][0]["Contacts"][0]["MaxMeshSizeEdge"] = edge_size

    def hp(c):
        c["Mesh"]["HPRefinement"] = {"Active": True, "Levels": 2, "Factor": 0.125}

    def best(c):
        # Finest static benchmark mesh. The MRI voxel size (0.5 mm) yields
        # ~18M elements, which overflows the .vtu writer (the file becomes
        # unreadable), so cap a bit coarser -- still by far the densest panel
        # while staying a readable file.
        c["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = max(2.0 * voxel_size, 1.0)

    strategies = [
        ("Results_VTA_default", default, 0, False),
        ("Results_VTA_fine", fine, 0, False),
        ("Results_VTA_edge_refinement", edge, 0, False),
        ("Results_VTA_material_refinement", default, 2, False),
        ("Results_VTA_hp_refinement", hp, 0, True),
        ("Results_VTA_best", best, 1, False),
    ]

    for out, mutate, material_steps, use_hp in strategies:
        c = json.loads(json.dumps(cfg))  # deep copy of a clean base config
        mutate(c)
        c["OutputPath"] = out
        _logger.info(f"=== generating {out} ===")
        export_strategy(c, mri_image, material_steps=material_steps, hp=use_hp)


if __name__ == "__main__":
    main()
