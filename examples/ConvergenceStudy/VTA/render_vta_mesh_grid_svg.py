"""3x2 grid of VTA mesh-refinement panels with an editable-text SVG.

Companion to the MeshRefinement figure
(``FiguresReproducibility/MeshRefinement/render_mesh_refinement_combined_svg.py``)
but for the VTA study meshes. Two differences on purpose:

1. **The mesh plots stay independent.** Every panel is rendered on its own
   (one PyVista screenshot per panel) and embedded as a separate ``<image>``
   in the SVG, so each mesh can be moved / swapped / re-scaled on its own in
   Inkscape/Illustrator. (The MeshRefinement script bakes all panels into one
   raster.) Only the text is emitted as real ``<text>`` (``svg.fonttype='none'``).

2. **A "Benchmark" panel.** The top-left panel shows the gold-standard
   ``Results_VTA_best`` mesh; the other five show the refinement strategies we
   already had. In the VTA study the Medtronic3389 lead is *tilted* (see
   ``base_settings.json`` -> Direction/TipPosition), so every panel is
   transformed into the electrode's own frame first, which makes the lead
   appear axis-aligned (vertical) exactly like the MeshRefinement panels.

This script NEVER generates meshes. It expects the ``Results_VTA_*`` result
dirs (each with a ``material.vtu``) to be precomputed by
``Medtronic3389/<scenario>/run_convergence_study.py`` and errors out listing
anything that is missing.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import contextlib
import io

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import ListedColormap
from matplotlib.patches import ConnectionPatch, Rectangle

pv.OFF_SCREEN = True
# keep text as editable <text> elements (system fonts), not outlined paths
plt.rcParams["svg.fonttype"] = "none"

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)

# --- what to render --------------------------------------------------------
ELECTRODE_NAME = "Medtronic3389"
# scenario subdir holding the precomputed Results_VTA_* dirs
SCENARIO = os.path.join("Medtronic3389", "case_grounding")
BASE_SETTINGS = "base_settings.json"  # for the (tilted) VTA lead placement
LEAD_VTM = "medtronic3389_lead_axis.vtm"  # axis-aligned lead cache

# (result dir, label) in grid (row-major). The first is the gold standard.
PANELS = [
    ("Results_VTA_best", "Benchmark"),
    ("Results_VTA_default", "Default (I)"),
    ("Results_VTA_fine", "Global fine (II-III)"),
    ("Results_VTA_edge_refinement", "Edge refinement (IV-VII)"),
    ("Results_VTA_material_refinement", "Material refinement (VIII-XII)"),
    ("Results_VTA_hp_refinement", "hp refinement (XIII-XIV)"),
]
BENCHMARK_INDEX = 0
GRID = (3, 2)  # rows, cols

# --- look (shared with the MeshRefinement figure) --------------------------
# material_real: 0 unknown, 1 gray matter, 2 white matter, 3 CSF, 4 blood
MATERIAL_COLORS = [
    "#dcdcdc",  # unknown (rare; sits right at the electrode)
    "#8b8f94",  # gray matter -> grey
    "#f4f5f7",  # white matter -> white
    "#2f6fb6",  # CSF -> blue
    "#b0413e",  # blood / encapsulation -> red (unused when thickness 0)
]
CMAP = ListedColormap(MATERIAL_COLORS)
CLIM = (0, 4)
EDGE_COLOR = "#2b2b2b"
LINE_WIDTH = 0.4
BODY_COLOR = (0.80, 0.80, 0.82)  # light grey lead body
CONTACT_COLOR = (0.22, 0.22, 0.24)  # dark grey contacts
# thin extracted-cell slab + soft lighting => 3-D shaded cells
SLAB_HALF = 0.1  # mm half-thickness
LIGHT_KW = {"lighting": True, "ambient": 0.6, "diffuse": 0.4, "specular": 0.0}

ORANGE = "#e8820c"  # benchmark highlight border
PAGE_BG = "#0a0a0a"  # black page, like the mock-up (gaps show through)

# --- camera (electrode frame: tip at origin, lead axis = +z) ---------------
FOCAL_OFFSET = 5.017  # mm above the tip along the lead axis (centres contacts)
PARALLEL_SCALE = 13.0  # half-height of the view in mm (zoomed out to show CSF)
FOCAL = (0.0, 0.0, FOCAL_OFFSET)

# --- output geometry (px) --------------------------------------------------
PANEL_W, PANEL_H = 560, 620
GAP = 12
MARGIN = 12
OUTFILE = "vta_mesh_grid"

# --- hp contact-edge magnifier (independent inset on the hp panel) ---------
# hp refinement is only visible right at the contact edges (singular corners),
# so the hp panel gets a zoomed contact-edge inset -- its own <image>.
HP_PANEL_INDEX = 5
HP_INSET_RESULT = "Results_VTA_hp_refinement"
HP_INSET_WINDOW = (300, 320)  # px
HP_INSET_SCALE = 0.4  # mm half-height of the zoomed contact-corner view


def electrode_frame():
    """Return (R, t) mapping VTA world coords into the electrode frame.

    In the electrode frame the lead tip is at the origin and the lead axis is
    +z, so a mesh transformed by ``R @ p + t`` renders as if the (tilted) VTA
    electrode were vertical -- matching the axis-aligned MeshRefinement panels.
    The slice plane is the electrode axis crossed with world +x.
    """
    with open(BASE_SETTINGS) as fp:
        e = json.load(fp)["Electrodes"][0]
    tip = e["TipPosition"]
    p0 = np.array([tip["x[mm]"], tip["y[mm]"], tip["z[mm]"]])
    d = np.array(
        [e["Direction"]["x[mm]"], e["Direction"]["y[mm]"], e["Direction"]["z[mm]"]]
    )
    ez = d / np.linalg.norm(d)
    # in-plane axis: world +x projected perpendicular to the lead axis
    ex = np.array([1.0, 0.0, 0.0]) - np.dot([1.0, 0.0, 0.0], ez) * ez
    if np.linalg.norm(ex) < 1e-6:  # lead ~parallel to x -> fall back to world +y
        ex = np.array([0.0, 1.0, 0.0]) - np.dot([0.0, 1.0, 0.0], ez) * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, ez)  # re-orthonormalise
    R = np.vstack([ex, ey, ez])
    return R, -R @ p0


def transform_matrix():
    """4x4 world->electrode-frame transform for pyvista.transform()."""
    R, t = electrode_frame()
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = t
    return M


def build_lead_parts():
    """Tessellate an axis-aligned Medtronic3389 by boundary; cache as VTM.

    Built with tip at the origin and axis = +z, so it lines up with the
    electrode-frame meshes. One PolyData block per boundary ('Body',
    'Contact_1', ...). Loads from the cache if present.
    """
    if os.path.exists(LEAD_VTM):
        return pv.read(LEAD_VTM)

    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, TaskManager

    import ossdbs

    settings = {
        "Mesh": {},
        "Electrodes": [
            {
                "Name": ELECTRODE_NAME,
                "TipPosition": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
                "Direction": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0},
                "Rotation[Degrees]": 0.0,
                "EncapsulationLayer": {"Thickness[mm]": 0.0},
            }
        ],
    }
    electrode = ossdbs.api.generate_electrodes(settings)[0]

    occgeo = OCCGeometry(electrode.geometry)
    with (  # netgen prints cosmetic warnings during surface meshing
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


def slab(result_dir, M):
    """Extract a thin slab of whole cells around the electrode plane.

    The mesh is first mapped into the electrode frame (``M``); two crinkle
    clips at y = +/- SLAB_HALF keep the tetrahedra whole (shaded 3-D facets)
    instead of cutting them flat.
    """
    mesh = pv.read(os.path.join(SCENARIO, result_dir, "material.vtu"))
    mesh = mesh.transform(M, inplace=False)
    return mesh.clip(
        normal=(0, 1, 0), origin=(0, -SLAB_HALF, 0), invert=False, crinkle=True
    ).clip(normal=(0, 1, 0), origin=(0, SLAB_HALF, 0), invert=True, crinkle=True)


def render_panel(result_dir, lead, M, path):
    """Render a single panel (no text) to its own PNG."""
    p = pv.Plotter(off_screen=True, window_size=[PANEL_W, PANEL_H], border=False)
    p.add_mesh(
        slab(result_dir, M),
        scalars="material_real",
        cmap=CMAP,
        clim=CLIM,
        show_edges=True,
        edge_color=EDGE_COLOR,
        line_width=LINE_WIDTH,
        show_scalar_bar=False,
        **LIGHT_KW,
    )
    add_lead(p, lead)
    p.set_background("white")
    p.enable_parallel_projection()
    p.camera_position = [
        (FOCAL[0], FOCAL[1] - 100.0, FOCAL[2]),
        FOCAL,
        (0.0, 0.0, 1.0),
    ]
    p.camera.parallel_scale = PARALLEL_SCALE
    p.screenshot(path)


def _frame_to_panel_px(x, z):
    """Pixel position within a panel of an electrode-frame (x, z) point."""
    half_h = PARALLEL_SCALE
    half_w = half_h * PANEL_W / PANEL_H
    px = (x - (FOCAL[0] - half_w)) / (2 * half_w) * PANEL_W
    py = ((FOCAL[2] + half_h) - z) / (2 * half_h) * PANEL_H
    return px, py


def render_hp_inset(M, focal, marker_ends, path):
    """Zoomed contact-corner view of the hp mesh, with an orange edge marker.

    No lead overlay (it would occlude the graded elements): the electrode side
    reads as empty white, the hp-graded tetrahedra hug the contact corner, and
    the orange marker traces the contact rim.
    """
    p = pv.Plotter(off_screen=True, window_size=list(HP_INSET_WINDOW), border=False)
    p.add_mesh(
        slab(HP_INSET_RESULT, M),
        scalars="material_real",
        cmap=CMAP,
        clim=CLIM,
        show_edges=True,
        edge_color="#1a1a1a",
        line_width=1.2,
        show_scalar_bar=False,
        **LIGHT_KW,
    )
    p.add_mesh(pv.Line(*marker_ends), color=ORANGE, line_width=6)
    p.set_background("white")
    p.enable_parallel_projection()
    p.camera_position = [(focal[0], focal[1] - 100.0, focal[2]), focal, (0.0, 0.0, 1.0)]
    p.camera.parallel_scale = HP_INSET_SCALE
    p.screenshot(path)


def _compose_hp_inset(fig, axes, lead, M, fig_w, fig_h, cols, rows):
    """Render the hp contact-edge magnifier and place it as its own image.

    Picks the contact nearest the view centre, renders a zoomed contact-edge
    view of the hp mesh, and composites it on the hp panel: a source box round
    the contact, orange connectors, and an independent inset <image>.
    """
    contacts = [k for k in lead.keys() if k.startswith("Contact")]
    z_centre = {k: 0.5 * (lead[k].bounds[4] + lead[k].bounds[5]) for k in contacts}
    target = min(contacts, key=lambda k: abs(z_centre[k] - FOCAL_OFFSET))
    zmin, zmax = lead[target].bounds[4], lead[target].bounds[5]
    lead_r = lead[target].bounds[1]  # xmax ~ lead radius in the x-z slice
    focal = (lead_r, 0.0, zmax)  # top rim corner, where hp grades most
    marker = ((lead_r + 0.02, 0.0, zmin), (lead_r + 0.02, 0.0, zmax))
    render_hp_inset(M, focal, marker, "_hp_inset.png")

    ax_hp = axes[HP_PANEL_INDEX]
    # source box on the hp panel matching the inset's mm extent
    half_w = HP_INSET_SCALE * HP_INSET_WINDOW[0] / HP_INSET_WINDOW[1]
    bx0, btop = _frame_to_panel_px(lead_r - half_w, focal[2] + HP_INSET_SCALE)
    bx1, bbot = _frame_to_panel_px(lead_r + half_w, focal[2] - HP_INSET_SCALE)
    ax_hp.add_patch(
        Rectangle(
            (bx0, btop),
            bx1 - bx0,
            bbot - btop,
            fill=False,
            edgecolor=ORANGE,
            linewidth=1.5,
            clip_on=False,
        )
    )

    # independent inset axes, anchored at the hp panel's bottom-left
    iw, ih = HP_INSET_WINDOW
    r_hp, c_hp = divmod(HP_PANEL_INDEX, cols)
    x0_hp = MARGIN + c_hp * (PANEL_W + GAP)
    y0_hp = MARGIN + (rows - 1 - r_hp) * (PANEL_H + GAP)
    pad = 12
    ax_in = fig.add_axes(
        ((x0_hp + pad) / fig_w, (y0_hp + pad) / fig_h, iw / fig_w, ih / fig_h)
    )
    ax_in.imshow(plt.imread("_hp_inset.png"))
    ax_in.set_xlim(0, iw)
    ax_in.set_ylim(ih, 0)
    ax_in.axis("off")
    ax_in.add_patch(
        Rectangle(
            (0, 0),
            iw - 1,
            ih - 1,
            fill=False,
            edgecolor=ORANGE,
            linewidth=3,
            clip_on=False,
        )
    )
    ax_in.text(8, 16, "hp @ contact edge", fontsize=11, color="black", va="top")

    # connectors: source-box bottom corners -> inset top corners
    for (box_x, box_y), (in_x, in_y) in [((bx0, bbot), (0, 0)), ((bx1, bbot), (iw, 0))]:
        fig.add_artist(
            ConnectionPatch(
                xyA=(in_x, in_y),
                coordsA=ax_in.transData,
                xyB=(box_x, box_y),
                coordsB=ax_hp.transData,
                color=ORANGE,
                linewidth=1.2,
            )
        )


def main():
    """Render every panel independently, then compose the grid SVG."""
    missing = [
        rd
        for rd, _ in PANELS
        if not os.path.exists(os.path.join(SCENARIO, rd, "material.vtu"))
    ]
    if missing:
        raise SystemExit(
            "Missing precomputed meshes (run run_convergence_study.py in "
            f"{SCENARIO} first); not recomputing here: {missing}"
        )

    M = transform_matrix()
    lead = build_lead_parts()

    tmp_pngs = []
    for idx, (result_dir, _label) in enumerate(PANELS):
        png = f"_panel_{idx}.png"
        render_panel(result_dir, lead, M, png)
        tmp_pngs.append(png)

    # --- compose: one independent <image> per panel + editable <text> ------
    rows, cols = GRID
    fig_w = 2 * MARGIN + cols * PANEL_W + (cols - 1) * GAP
    fig_h = 2 * MARGIN + rows * PANEL_H + (rows - 1) * GAP
    dpi = 100
    fig = plt.figure(figsize=(fig_w / dpi, fig_h / dpi), dpi=dpi, facecolor=PAGE_BG)

    axes = []
    for idx, (_result_dir, label) in enumerate(PANELS):
        r, c = divmod(idx, cols)
        x0 = MARGIN + c * (PANEL_W + GAP)
        # matplotlib figure fraction has y=0 at the bottom -> flip the row
        y0 = MARGIN + (rows - 1 - r) * (PANEL_H + GAP)
        ax = fig.add_axes((x0 / fig_w, y0 / fig_h, PANEL_W / fig_w, PANEL_H / fig_h))
        ax.imshow(plt.imread(tmp_pngs[idx]))
        ax.set_xlim(0, PANEL_W)
        ax.set_ylim(PANEL_H, 0)
        ax.axis("off")
        axes.append(ax)

        is_bench = idx == BENCHMARK_INDEX
        ax.add_patch(
            Rectangle(
                (0, 0),
                PANEL_W,
                PANEL_H,
                fill=False,
                edgecolor=ORANGE if is_bench else "#000000",
                linewidth=6 if is_bench else 2,
                clip_on=False,
            )
        )
        ax.text(
            18,
            30,
            label,
            fontsize=24 if is_bench else 18,
            color="black",
            va="top",
            ha="left",
            weight="bold",
        )

    _compose_hp_inset(fig, axes, lead, M, fig_w, fig_h, cols, rows)

    fig.savefig(f"{OUTFILE}.svg", facecolor=PAGE_BG)
    fig.savefig(f"{OUTFILE}_svg_preview.png", dpi=dpi, facecolor=PAGE_BG)
    plt.close(fig)

    for png in [*tmp_pngs, "_hp_inset.png"]:
        os.remove(png)
    print(f"wrote {OUTFILE}.svg and {OUTFILE}_svg_preview.png")


if __name__ == "__main__":
    main()
