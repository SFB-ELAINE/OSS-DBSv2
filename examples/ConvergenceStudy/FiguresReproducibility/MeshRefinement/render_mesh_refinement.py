"""Render a mesh-refinement comparison figure with PyVista.

For each representative refinement strategy this slices the FEM mesh
(``Results_PAM_<strategy>/material.vtu``) through the electrode axis, shows
the cut triangulation (mesh wireframe) over the tissue, and overlays the
electrode (light-grey body, dark-grey contacts). The panels are assembled
into one single-row overview that highlights how the strategies refine the
mesh near the electrode.

Run ``generate_meshes.py`` first to produce the ``material.vtu`` files.
"""

import contextlib
import io
import json
import os

import numpy as np
import pyvista as pv
from matplotlib.colors import ListedColormap

pv.OFF_SCREEN = True

CONFIG = "oss_dbs_parameters.json"
ELECTRODE_NAME = "Medtronic3389"
LEAD_VTM = "medtronic3389_lead.vtm"  # per-boundary surface mesh cache

# --- electrode / slice geometry (Medtronic3389, axis = +z) ---
ELECTRODE_X = 12.956271802141353
ELECTRODE_Y = -9.901870551007098
# slice plane through the lead axis -> shows the x-z cross-section
SLICE_ORIGIN = (ELECTRODE_X, ELECTRODE_Y, 0.0)
SLICE_NORMAL = (0.0, 1.0, 0.0)
# camera: look along -y onto the x-z plane, zoomed to the contact region
FOCAL = (ELECTRODE_X, ELECTRODE_Y, -7.0)
PARALLEL_SCALE = 9.5  # half-height of the view in mm

# material_real: 0 unknown, 1 gray matter, 2 white matter, 3 CSF
MATERIAL_COLORS = [
    "#dcdcdc",  # unknown (rare; sits right at the electrode)
    "#8b8f94",  # gray matter -> grey
    "#f4f5f7",  # white matter -> white
    "#2f6fb6",  # CSF -> blue
]
CMAP = ListedColormap(MATERIAL_COLORS)
CLIM = (0, 3)

EDGE_COLOR = "#2b2b2b"  # dark, so the wireframe shows on white matter too
LINE_WIDTH = 0.4

BODY_COLOR = (0.80, 0.80, 0.82)  # light grey lead body
CONTACT_COLOR = (0.22, 0.22, 0.24)  # dark grey contacts

# (result dir, panel label) — canonical numerals; default_meshsize is the
# neuron mesh-size strategy (excluded from the paper convergence figures).
PANELS = [
    ("Results_PAM_default", "Default (I)"),
    ("Results_PAM_fine_edge_refinement", "Edge refinement (V)"),
    ("Results_PAM_material_refinement", "Material refinement (VIII)"),
    ("Results_PAM_default_meshsize", "Pathway mesh size (neuron)"),
]
WINDOW = (3200, 850)
OUTFILE = "mesh_refinement_overview"


def build_lead_parts():
    """Tessellate the electrode by named boundary; cache as a MultiBlock VTM.

    Returns a MultiBlock with one PolyData block per boundary ('Body',
    'Contact_1', ...). Loads from the cache if present.
    """
    if os.path.exists(LEAD_VTM):
        return pv.read(LEAD_VTM)

    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, TaskManager

    import ossdbs

    with open(CONFIG) as fp:
        settings = json.load(fp)
    settings["Electrodes"][0]["Name"] = ELECTRODE_NAME
    electrode = ossdbs.api.generate_electrodes(settings)[0]

    occgeo = OCCGeometry(electrode.geometry)
    # netgen prints cosmetic warnings during surface meshing
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
        TaskManager(),
    ):
        ngmesh = occgeo.GenerateMesh(maxh=0.3)

    m = Mesh(ngmesh)
    bnds = m.GetBoundaries()  # names, 0-indexed; element .index is 1-based
    pts = np.array([p.p for p in ngmesh.Points()])

    tris_by_name = {}
    for el in ngmesh.Elements2D():
        bname = bnds[el.index - 1]
        tris_by_name.setdefault(bname, []).append([v.nr - 1 for v in el.vertices])

    remap = np.empty(len(pts), dtype=np.int64)
    mb = pv.MultiBlock()
    for name in sorted(tris_by_name):
        tris = np.array(tris_by_name[name], dtype=np.int64)
        used = np.unique(tris)
        remap[used] = np.arange(len(used), dtype=np.int64)
        faces = np.hstack(
            [np.full((len(tris), 1), 3, dtype=np.int64), remap[tris]]
        ).ravel()
        mb.append(pv.PolyData(pts[used], faces), name=name)

    mb.save(LEAD_VTM)
    print(f"built {LEAD_VTM} with boundaries: {sorted(tris_by_name)}")
    return mb


def add_lead(plotter, lead):
    """Overlay the lead: light-grey body, dark-grey contacts."""
    for name in lead.keys():
        color = CONTACT_COLOR if name.startswith("Contact") else BODY_COLOR
        plotter.add_mesh(
            lead[name],
            color=color,
            smooth_shading=True,
            split_sharp_edges=True,
            show_edges=False,
        )


def slice_mesh(result_dir):
    """Slice a strategy's material mesh through the electrode axis.

    ``generate_triangles=False`` keeps each tetrahedron's natural cut
    polygon (a triangle or quad) instead of triangulating it, which avoids
    the spurious diagonal sliver edges a plain slice introduces.
    """
    mesh = pv.read(os.path.join(result_dir, "material.vtu"))
    return mesh.slice(
        normal=SLICE_NORMAL, origin=SLICE_ORIGIN, generate_triangles=False
    )


def main():
    """Assemble the single-row mesh-refinement overview."""
    missing = [d for d, _ in PANELS if not os.path.exists(f"{d}/material.vtu")]
    if missing:
        raise SystemExit(f"Missing meshes (run generate_meshes.py): {missing}")

    lead = build_lead_parts()
    p = pv.Plotter(
        off_screen=True, shape=(1, len(PANELS)), window_size=list(WINDOW), border=False
    )
    for idx, (result_dir, label) in enumerate(PANELS):
        p.subplot(0, idx)
        sl = slice_mesh(result_dir)
        p.add_mesh(
            sl,
            scalars="material_real",
            cmap=CMAP,
            clim=CLIM,
            show_edges=True,
            edge_color=EDGE_COLOR,
            line_width=LINE_WIDTH,
            lighting=False,
            show_scalar_bar=False,
        )
        add_lead(p, lead)
        p.add_text(
            label, font_size=18, color="black", position="upper_left", shadow=True
        )
        p.set_background("white")
        p.enable_parallel_projection()
        p.camera_position = [
            (FOCAL[0], FOCAL[1] - 100.0, FOCAL[2]),
            FOCAL,
            (0.0, 0.0, 1.0),
        ]
        p.camera.parallel_scale = PARALLEL_SCALE

    p.screenshot(f"{OUTFILE}.png")
    print(f"wrote {OUTFILE}.png")


if __name__ == "__main__":
    main()
