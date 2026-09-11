"""Generate the Benchmark panel's mesh at the true voxel (0.5 mm) resolution.

The finest static benchmark mesh (voxel-size limit + material refinement) has
~18M elements; writing the whole mesh to ``.vtu`` overflows the ASCII writer
(the >2 GB file comes back corrupt). But the figure only shows a thin slice, so
this builds the full-resolution mesh and keeps only the render slab *before*
writing -- a small, valid binary ``.vtu`` with the real 0.5 mm density in the
cross-section that is actually drawn.

No FEM solve: material labels come from an L2(order 0) projection of the
material distribution (local, one constant per element), matching what
``export_material_distribution_to_vtk`` samples. Run from this directory
(Medtronic3389/case_grounding); reuses ``generate_vta_meshes`` for the setup.

Writes ``Results_VTA_best/material.vtu`` (overwriting the coarse placeholder).
"""

import json
import logging
import os

import generate_vta_meshes as g
import ngsolve
import numpy as np
import pyvista as pv

import ossdbs

ossdbs.set_logger(logging.INFO)
_logger = logging.getLogger("ossdbs")

OUT_DIR = "Results_VTA_best"
# world-space half-thickness of the exported slab (mm). Must exceed the
# render's SLAB_HALF (0.1) so its crinkle clip has whole cells to keep.
SLAB_HALF = 0.3


def slab_plane():
    """(origin p0, normal e_y) of the render slice plane, in world coords.

    Mirrors ``render_vta_mesh_grid_svg.electrode_frame``: the slice plane
    contains the lead axis and is spanned so its normal is the electrode
    frame's +y.
    """
    with open("../../base_settings.json") as fp:
        e = json.load(fp)["Electrodes"][0]
    tip = e["TipPosition"]
    p0 = np.array([tip["x[mm]"], tip["y[mm]"], tip["z[mm]"]])
    d = np.array(
        [e["Direction"]["x[mm]"], e["Direction"]["y[mm]"], e["Direction"]["z[mm]"]]
    )
    ez = d / np.linalg.norm(d)
    ex = np.array([1.0, 0.0, 0.0]) - np.dot([1.0, 0.0, 0.0], ez) * ez
    if np.linalg.norm(ex) < 1e-6:
        ex = np.array([0.0, 1.0, 0.0]) - np.dot([0.0, 1.0, 0.0], ez) * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    ey /= np.linalg.norm(ey)
    return p0, ey


def main():
    """Mesh at voxel resolution, keep the slab, write a small binary .vtu."""
    p0, ey = slab_plane()

    cfg = g.base_config()
    mri_image, _ = ossdbs.load_images(cfg)
    voxel_size = float(min(mri_image.voxel_sizes))
    cfg["Mesh"]["MeshingHypothesis"]["MaxMeshSize"] = voxel_size  # true 0.5 mm
    cfg["OutputPath"] = OUT_DIR

    geometry, brain_model = g.build_geometry(cfg, mri_image)
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
        vc.refine_mesh_by_material(1)  # match the study best (1 material step)
        mesh = vc.mesh.ngsolvemesh
        _logger.info(f"{OUT_DIR}: {mesh.ne} elements (full voxel mesh)")

        # per-element material label (local L2 projection -> no solve)
        material_cf = vc.conductivity_cf.material_distribution(vc.mesh)
        gf = ngsolve.GridFunction(ngsolve.L2(mesh, order=0))
        gf.Set(material_cf)
        material = np.asarray(gf.vec.FV()).real

    ngm = mesh.ngmesh
    pts = np.asarray(ngm.Coordinates())
    conn = np.asarray(ngm.Elements3D().NumPy()["nodes"])[:, :4].astype(np.int64) - 1

    # keep tets whose signed-distance range overlaps the slab [-h, h]
    sd = (pts - p0) @ ey
    sdc = sd[conn]
    keep = (sdc.min(axis=1) <= SLAB_HALF) & (sdc.max(axis=1) >= -SLAB_HALF)
    conn = conn[keep]
    material = material[keep]
    _logger.info(f"slab: {len(conn)} of {len(sd)} points-worth cells kept")

    used = np.unique(conn)
    remap = np.full(len(pts), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    cells = np.hstack([np.full((len(conn), 1), 4, dtype=np.int64), remap[conn]]).ravel()
    celltypes = np.full(len(conn), pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, celltypes, pts[used])
    grid.cell_data["material_real"] = material

    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, "material.vtu")
    grid.save(out)  # pyvista writes a valid binary .vtu
    _logger.info(f"wrote {out}: {grid.n_cells} cells, {grid.n_points} points")


if __name__ == "__main__":
    main()
