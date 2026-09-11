"""Generate the Benchmark mesh: the exact study-best (voxel-0.5 mm + material
+ adaptive mesh refinement).

Builds the study-best base (voxel-0.5 mm + 1 material step) and runs
``run_full_analysis`` with adaptive mesh refinement (AMR) for 3 steps, with no
VTA point analysis and no VTK solution export. Logs the element count before
and after each refinement, then extracts the render slab of the final mesh to
``Results_VTA_best/material.vtu`` (this is the figure's Benchmark panel).

Empirically, AMR adds only ~0.1% elements here (the 0.5 mm base is already
converged), so the benchmark is essentially uniformly ultra-fine.

Run from this directory (Medtronic3389/case_grounding) after
``generate_vta_meshes.py``.
"""

import logging
import os

import generate_best_slab as gbs
import generate_vta_meshes as g
import ngsolve
import numpy as np
import pyvista as pv

import ossdbs

ossdbs.set_logger(logging.INFO)
_logger = logging.getLogger("ossdbs")

OUT_DIR = "Results_VTA_best"
AMR_STEPS = 3


def main():
    """Build the study-best mesh with AMR and write its render slab."""
    p0, ey = gbs.slab_plane()

    cfg = g.base_config()
    mri_image, _ = ossdbs.load_images(cfg)
    voxel_size = float(min(mri_image.voxel_sizes))
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = voxel_size  # true 0.5 mm
    cfg["OutputPath"] = OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)

    geometry, brain_model = g.build_geometry(cfg, mri_image)
    solver = ossdbs.api.prepare_solver(cfg)
    dielectric_properties = ossdbs.api.prepare_dielectric_properties(cfg)
    conductivity = ossdbs.ConductivityCF(
        mri_image,
        brain_model.brain_region,
        dielectric_properties,
        cfg["MaterialDistribution"]["MRIMapping"],
        geometry.encapsulation_layers,
        complex_data=cfg["EQSMode"],
    )
    signal = ossdbs.api.prepare_stimulation_signal(cfg)

    with ngsolve.TaskManager():
        vc = ossdbs.api.prepare_volume_conductor_model(
            cfg, geometry, conductivity, solver
        )
        vc.refine_mesh_by_material(1)  # study-best base: 1 material step
        _logger.info(f"best base elements: {vc.mesh.ngsolvemesh.ne}")

        # force AMR_STEPS refinements (tolerance 0 -> never satisfied)
        vc.run_full_analysis(
            frequency_domain_signal=signal,
            compute_impedance=False,
            export_vtk=False,
            point_models=[],
            adaptive_mesh_refinement_settings={
                "Active": True,
                "MaxIterations": AMR_STEPS,
                "ErrorTolerance": 0.0,
            },
        )
        mesh = vc.mesh.ngsolvemesh
        _logger.info(f"best final elements: {mesh.ne}")

        material_cf = vc.conductivity_cf.material_distribution(vc.mesh)
        gf = ngsolve.GridFunction(ngsolve.L2(mesh, order=0))
        gf.Set(material_cf)
        material = np.asarray(gf.vec.FV()).real

    ngm = mesh.ngmesh
    pts = np.asarray(ngm.Coordinates())
    conn = np.asarray(ngm.Elements3D().NumPy()["nodes"])[:, :4].astype(np.int64) - 1
    sd = (pts - p0) @ ey
    sdc = sd[conn]
    keep = (sdc.min(axis=1) <= gbs.SLAB_HALF) & (sdc.max(axis=1) >= -gbs.SLAB_HALF)
    conn = conn[keep]
    material = material[keep]
    used = np.unique(conn)
    remap = np.full(len(pts), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    cells = np.hstack([np.full((len(conn), 1), 4, dtype=np.int64), remap[conn]]).ravel()
    celltypes = np.full(len(conn), pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, celltypes, pts[used])
    grid.cell_data["material_real"] = material
    out = os.path.join(OUT_DIR, "material.vtu")
    grid.save(out)
    _logger.info(f"wrote {out}: {grid.n_cells} cells, {grid.n_points} points")


if __name__ == "__main__":
    main()
