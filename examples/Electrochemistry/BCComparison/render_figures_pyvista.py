"""PyVista counterpart of render_figures.py.

Same six runs, same slice plane, same camera, same fixed scalar
ranges — just rendered with VTK + PyVista so it works in the venv
without a ParaView install. Run as:

    python render_figures_pyvista.py                       # all six
    python render_figures_pyvista.py A1_floating
    python render_figures_pyvista.py A1_floating B1_dirichlet

Outputs sit alongside the pvbatch ones, with `_pv` suffix for direct
side-by-side comparison:

    Figures_<run_id>/potential_slice_pv.png
    Figures_<run_id>/E_field_slice_pv.png

Lead overlay
------------
The Vercise lead is rendered with two face groups (Body = silver,
Contact_1…8 = gold) built from the ngsolve surface tessellation, plus
pyvista-cad B-rep edges from the STEP export. Both artefacts are cached
on disk so subsequent runs are fast.  Edge appearance is controlled by
EDGE_LINE_WIDTH and EDGE_AS_TUBES at the top of the file.
"""

import contextlib
import io
import json
import os
import string
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv
import pyvista_cad  # noqa: F401 — registers pv.read STEP support + plotter.cad
from matplotlib.colors import ListedColormap

# Slice plane + camera are derived from base_settings.json by
# generate_paraview_camera.py. Regenerate after editing the JSON.
from paraview_camera import (
    CAMERA_FOCAL_POINT,
    CAMERA_PARALLEL_SCALE,
    CAMERA_POSITION,
    CAMERA_VIEW_UP,
    SLICE_NORMAL,
    SLICE_ORIGIN,
)

ALL_RUNS = [
    "A1_floating",
    "A2_insulating",
    "A3_mild_Z",
    "B1_dirichlet",
    "B2_mild_Z",
    "B3_strong_Z",
]

# ---- Lead geometry files (auto-generated, safe to delete to rebuild) ----
LEAD_STEP = "vercise_lead.stp"  # STEP export for B-rep edges + external use
LEAD_VTM = "vercise_lead.vtm"  # per-boundary surface mesh cache

# ---- Lead appearance ----
BODY_COLOR = (0.78, 0.78, 0.80)  # titanium / silver
CONTACT_COLOR = (0.90, 0.78, 0.20)  # gold / platinum  (non-active fallback)
LEAD_OPACITY = 0.85
EDGE_COLOR = "dimgray"

# Per-contact overrides: red = higher voltage, blue = lower voltage.
# Contacts not listed here fall back to CONTACT_COLOR.
ACTIVE_CONTACT_COLORS: dict[str, tuple[float, float, float]] = {
    "Contact_2": (0.18, 0.42, 0.80),  # 0 V — blue
    "Contact_4": (0.82, 0.22, 0.18),  # 1 V — red
}
EDGE_LINE_WIDTH = 3.0  # ← increase for thicker edges
EDGE_AS_TUBES = True  # ← True gives round-capped 3-D tubes, False = flat lines

# ---- Field data ----
POTENTIAL_RANGE = (0.0, 1.0)
E_FIELD_RANGE = (1.0, 300.0)

MATERIAL_LABELS = ["Unknown", "Gray matter", "White matter", "CSF"]
MATERIAL_COLORS = [
    (1.00, 1.00, 1.00),
    (0.85, 0.70, 0.70),
    (0.95, 0.92, 0.78),
    (0.70, 0.85, 0.95),
]
MATERIAL_RANGE = (-0.5, 3.5)

IMAGE_SIZE = (1200, 1200)
POTENTIAL_CMAP = "coolwarm"
E_FIELD_CMAP = "inferno"

# Scalar bar: vertical strip on the left side.  All call sites start from
# this base dict and add a "title" (and optionally "fmt") key.
_SBAR = {
    "vertical": True,
    "position_x": 0.03,  # left edge of the widget (normalised 0..1)
    "position_y": 0.25,  # bottom edge of the widget
    "width": 0.10,  # total widget width (colour strip + labels)
    "height": 0.50,  # total widget height
    "n_labels": 5,
    "color": "white",
    "title_font_size": 26,
    "label_font_size": 22,
}

pv.OFF_SCREEN = True


# ---------------------------------------------------------------------------
# Lead geometry helpers
# ---------------------------------------------------------------------------


def _export_lead_step() -> None:
    """Export the Vercise lead as STEP (once; useful for ParaView / FreeCAD)."""
    if os.path.exists(LEAD_STEP):
        return
    import ossdbs

    with open("base_settings.json") as f:
        settings = json.load(f)
    electrodes = ossdbs.generate_electrodes(settings)
    electrodes[0].geometry.WriteStep(LEAD_STEP)
    print(f"[ok]   lead: exported {LEAD_STEP}")


def _build_lead_parts() -> pv.MultiBlock | None:
    """Tessellate the electrode by named boundary and cache as VTM.

    Returns a MultiBlock with one PolyData block per boundary name:
    'Body', 'Contact_1', …, 'Contact_8'.  Loads from LEAD_VTM if it
    already exists.
    """
    if os.path.exists(LEAD_VTM):
        return pv.read(LEAD_VTM)

    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, TaskManager

    import ossdbs

    with open("base_settings.json") as f:
        settings = json.load(f)

    electrode = ossdbs.generate_electrodes(settings)[0]
    shape = electrode.geometry

    occgeo = OCCGeometry(shape)
    # Suppress netgen console output during meshing (cosmetic warnings only).
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        with TaskManager():
            ngmesh = occgeo.GenerateMesh(maxh=0.3)

    m = Mesh(ngmesh)
    bnds = m.GetBoundaries()  # tuple of names, 0-indexed; element .index is 1-based

    pts_arr = np.array([p.p for p in ngmesh.Points()])  # (N, 3)

    tris_by_name: dict[str, list] = {}
    for el in ngmesh.Elements2D():
        bname = bnds[el.index - 1]
        tris_by_name.setdefault(bname, []).append([v.nr - 1 for v in el.vertices])

    # Remap index helper: maps global vertex indices to a compact local array.
    remap_buf = np.empty(len(pts_arr), dtype=np.int64)

    mb = pv.MultiBlock()
    for name in sorted(tris_by_name):
        tris = np.array(tris_by_name[name], dtype=np.int64)
        used = np.unique(tris)
        remap_buf[used] = np.arange(len(used), dtype=np.int64)
        local_pts = pts_arr[used]
        local_tris = remap_buf[tris]
        faces_flat = np.hstack(
            [np.full((len(local_tris), 1), 3, dtype=np.int64), local_tris]
        ).ravel()
        mb.append(pv.PolyData(local_pts, faces_flat), name=name)

    mb.save(LEAD_VTM)
    print(
        f"[ok]   lead: built surface mesh {LEAD_VTM} "
        f"({len(tris_by_name)} boundary groups)"
    )
    return mb


def _load_brep_edges() -> pv.PolyData | None:
    """Return the 66 exact B-rep edge polylines from the STEP file."""
    if not os.path.exists(LEAD_STEP):
        return None
    poly = pv.read(LEAD_STEP)[0]  # single-block PolyData with cached B-rep
    cv = poly.cad.cad_view()
    edges: pv.PolyData = cv["edges"]  # type: ignore[assignment]
    return edges if edges.n_cells else None


def _add_lead(
    plotter: pv.Plotter,
    lead_parts: pv.MultiBlock | None,
    lead_edges: pv.PolyData | None,
) -> None:
    """Overlay the lead on *plotter*.

    Body faces are rendered in BODY_COLOR, each contact in CONTACT_COLOR.
    B-rep edges (the exact topological feature curves) are drawn on top.
    """
    if lead_parts is None:
        return

    for name in lead_parts.keys():
        block = lead_parts[name]
        if not name.startswith("Contact"):
            color = BODY_COLOR
        else:
            color = ACTIVE_CONTACT_COLORS.get(name, CONTACT_COLOR)
        plotter.add_mesh(
            block,
            color=color,
            opacity=LEAD_OPACITY,
            smooth_shading=True,
            split_sharp_edges=True,
            show_edges=False,
        )

    if lead_edges is not None and lead_edges.n_cells:
        plotter.add_mesh(
            lead_edges,
            color=EDGE_COLOR,
            line_width=EDGE_LINE_WIDTH,
            render_lines_as_tubes=EDGE_AS_TUBES,
            pickable=False,
            lighting=False,
        )


# ---------------------------------------------------------------------------
# Shared rendering helpers
# ---------------------------------------------------------------------------


def _setup_camera(plotter: pv.Plotter) -> None:
    plotter.enable_parallel_projection()
    plotter.camera_position = [
        list(CAMERA_POSITION),
        list(CAMERA_FOCAL_POINT),
        list(CAMERA_VIEW_UP),
    ]
    plotter.camera.parallel_scale = CAMERA_PARALLEL_SCALE
    plotter.background_color = "white"


def _slice(path: str) -> pv.PolyData:
    return pv.read(path).slice(normal=SLICE_NORMAL, origin=SLICE_ORIGIN)


def _save(plotter: pv.Plotter, path: str) -> None:
    plotter.show(screenshot=path, auto_close=False)
    plotter.save_graphic(path.replace(".png", ".svg"))
    plotter.close()


# ---------------------------------------------------------------------------
# Per-run renderers
# ---------------------------------------------------------------------------


def render_run(
    run_id: str,
    lead_parts: pv.MultiBlock | None = None,
    lead_edges: pv.PolyData | None = None,
    show_scalar_bar: bool = True,
) -> None:
    """Run the renderer."""
    result_dir = f"Results_{run_id}"
    out_dir = f"Figures_{run_id}"
    os.makedirs(out_dir, exist_ok=True)

    pot_path = os.path.join(result_dir, "potential.vtu")
    e_path = os.path.join(result_dir, "E-field.vtu")
    if not (os.path.exists(pot_path) and os.path.exists(e_path)):
        print(f"[skip] {run_id}: missing potential.vtu or E-field.vtu")
        return

    suffix = "_pv" if show_scalar_bar else "_pv_nobar"

    # ---- potential slice ----
    pot_slice = _slice(pot_path)
    p = pv.Plotter(off_screen=True, window_size=list(IMAGE_SIZE))
    p.add_mesh(
        pot_slice,
        scalars="potential_real",
        cmap=POTENTIAL_CMAP,
        clim=POTENTIAL_RANGE,
        show_scalar_bar=show_scalar_bar,
        scalar_bar_args={**_SBAR, "title": "Electric potential / V", "fmt": "%.2f"},
    )
    _add_lead(p, lead_parts, lead_edges)
    _setup_camera(p)
    _save(p, os.path.join(out_dir, f"potential_slice{suffix}.png"))

    # ---- E-field magnitude slice ----
    e_slice = _slice(e_path)
    e_slice["|E|"] = 1000.0 * np.linalg.norm(e_slice["E_field_real"], axis=1)
    p = pv.Plotter(off_screen=True, window_size=list(IMAGE_SIZE))
    p.add_mesh(
        e_slice,
        scalars="|E|",
        cmap=E_FIELD_CMAP,
        clim=E_FIELD_RANGE,
        log_scale=True,
        show_scalar_bar=show_scalar_bar,
        scalar_bar_args={**_SBAR, "title": "E / V·m⁻¹", "fmt": "%.0f"},
    )
    _add_lead(p, lead_parts, lead_edges)
    _setup_camera(p)
    _save(p, os.path.join(out_dir, f"E_field_slice{suffix}.png"))

    print(f"[ok]   {run_id}: wrote {out_dir}/*_slice{suffix}.{{png,svg}}")


def render_material(
    run_id: str = "A1_floating",
    lead_parts: pv.MultiBlock | None = None,
    lead_edges: pv.PolyData | None = None,
) -> None:
    """Render the material / tissue distribution."""
    mat_path = os.path.join(f"Results_{run_id}", "material.vtu")
    if not os.path.exists(mat_path):
        print(f"[skip] material: missing {mat_path}")
        return
    out_dir = "Figures_material"
    os.makedirs(out_dir, exist_ok=True)

    mat_slice = _slice(mat_path)
    cmap = ListedColormap(MATERIAL_COLORS)
    p = pv.Plotter(off_screen=True, window_size=list(IMAGE_SIZE))
    p.add_mesh(
        mat_slice,
        scalars="material_real",
        cmap=cmap,
        clim=MATERIAL_RANGE,
        scalar_bar_args={
            **_SBAR,
            "title": "Material",
            "n_labels": len(MATERIAL_LABELS),
            "fmt": "%.0f",
        },
    )
    _add_lead(p, lead_parts, lead_edges)
    _setup_camera(p)
    _save(p, os.path.join(out_dir, "material_slice_pv.png"))
    print(f"[ok]   material: wrote {out_dir}/material_slice_pv.{{png,svg}}")


# ---------------------------------------------------------------------------
# Grid assembly
# ---------------------------------------------------------------------------

# Row 0 = Panel A variants, row 1 = Panel B variants.
_GRID_RUNS = [
    ["A1_floating", "A2_insulating", "A3_mild_Z"],
    ["B1_dirichlet", "B2_mild_Z", "B3_strong_Z"],
]


_SVG_NS = "http://www.w3.org/2000/svg"
_XLINK_NS = "http://www.w3.org/1999/xlink"

# Grid layout: gap between panels (in source SVG units = 1200 pt each).
_GRID_GAP = 20  # px between cells
_CELL_SIZE = 1200  # px — matches IMAGE_SIZE


def _combine_svg_grid(
    figure_type: str,
    slug_with_bar: str,
    slug_no_bar: str,
    out_path: str,
) -> None:
    """Assemble a 2x3 SVG grid from per-run SVG files.

    Panel A (index 0) uses the with-colorbar SVG; all others use the
    no-colorbar variant.  The SVG text elements (tick labels, title) from
    pyvista are preserved as editable vector text nodes.

    Labels A-F are added as SVG ``<text>`` elements in the top-left corner
    of each panel, so they too remain fully editable.
    """
    # Register namespaces so ElementTree does not emit ns0: prefixes.
    ET.register_namespace("", _SVG_NS)
    ET.register_namespace("xlink", _XLINK_NS)

    n_cols, n_rows = 3, 2
    cell = _CELL_SIZE
    gap = _GRID_GAP
    total_w = n_cols * cell + (n_cols - 1) * gap
    total_h = n_rows * cell + (n_rows - 1) * gap

    # Build root element without explicit xmlns attrs — ET writes them
    # automatically based on the tag URI and registered namespaces.
    root = ET.Element(
        f"{{{_SVG_NS}}}svg",
        attrib={
            "version": "1.1",
            "width": f"{total_w}",
            "height": f"{total_h}",
            "viewBox": f"0 0 {total_w} {total_h}",
        },
    )

    for idx in range(n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        run_id = _GRID_RUNS[row][col]

        is_first = idx == 0
        slug = slug_with_bar if is_first else slug_no_bar
        svg_path = os.path.join(f"Figures_{run_id}", f"{slug}.svg")

        tx = col * (cell + gap)
        ty = row * (cell + gap)

        if os.path.exists(svg_path):
            tree = ET.parse(svg_path)
            src_root = tree.getroot()
            # Wrap the entire source SVG content in a translated <g>.
            g = ET.SubElement(
                root,
                f"{{{_SVG_NS}}}g",
                attrib={"transform": f"translate({tx},{ty})"},
            )
            # Copy all child elements from the source SVG's root.
            for child in src_root:
                g.append(child)
        else:
            # Placeholder rectangle when the SVG is missing.
            ET.SubElement(
                root,
                f"{{{_SVG_NS}}}rect",
                attrib={
                    "x": str(tx),
                    "y": str(ty),
                    "width": str(cell),
                    "height": str(cell),
                    "fill": "#e0e0e0",
                },
            )

        # Panel letter label (fully editable SVG text).
        label_x = tx + int(0.02 * cell)
        label_y = ty + int(0.06 * cell)
        # White backing rectangle for legibility.
        pad = 10
        ET.SubElement(
            root,
            f"{{{_SVG_NS}}}rect",
            attrib={
                "x": str(label_x - pad),
                "y": str(label_y - int(0.05 * cell) - pad),
                "width": str(int(0.065 * cell) + 2 * pad),
                "height": str(int(0.065 * cell) + 2 * pad),
                "fill": "white",
                "fill-opacity": "0.75",
            },
        )
        ET.SubElement(
            root,
            f"{{{_SVG_NS}}}text",
            attrib={
                "x": str(label_x),
                "y": str(label_y),
                "font-family": "sans-serif",
                "font-size": str(int(0.06 * cell)),
                "font-weight": "bold",
                "fill": "black",
            },
        ).text = string.ascii_uppercase[idx]

    tree_out = ET.ElementTree(root)
    ET.indent(tree_out, space="  ")
    tree_out.write(out_path, xml_declaration=True, encoding="unicode")
    print(f"[ok]   grid: wrote {out_path}")


def render_grid(
    figure_type: str,
    lead_parts: pv.MultiBlock | None = None,
    lead_edges: pv.PolyData | None = None,
) -> None:
    """Render no-colorbar variants for panels B-F, then assemble SVG grid.

    Panel A (A1_floating) keeps its colorbar; all other panels are
    re-rendered without one so the grid has a single shared colorbar.
    """
    slug_root = "potential_slice" if figure_type == "potential" else "E_field_slice"
    slug_bar = f"{slug_root}_pv"
    slug_nobar = f"{slug_root}_pv_nobar"

    for idx, row in enumerate(_GRID_RUNS):
        for jdx, run_id in enumerate(row):
            is_first = idx == 0 and jdx == 0
            if is_first:
                # already rendered with bar; skip
                continue
            nobar_path = os.path.join(f"Figures_{run_id}", f"{slug_nobar}.svg")
            if not os.path.exists(nobar_path):
                render_run(run_id, lead_parts, lead_edges, show_scalar_bar=False)

    _combine_svg_grid(
        figure_type,
        slug_bar,
        slug_nobar,
        f"{figure_type}_grid.svg",
    )


def render_grids(
    lead_parts: pv.MultiBlock | None = None,
    lead_edges: pv.PolyData | None = None,
) -> None:
    """Render SVG grids (needed for paper figure)."""
    render_grid("potential", lead_parts, lead_edges)
    render_grid("E_field", lead_parts, lead_edges)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str]) -> None:
    """Main function."""
    _export_lead_step()  # write STEP once (for external tools)
    lead_parts = _build_lead_parts()
    lead_edges = _load_brep_edges()

    runs = argv[1:] if len(argv) > 1 else ALL_RUNS
    for run_id in runs:
        render_run(run_id, lead_parts, lead_edges)
    if len(argv) <= 1:
        render_material(lead_parts=lead_parts, lead_edges=lead_edges)
        render_grids(lead_parts=lead_parts, lead_edges=lead_edges)


if __name__ == "__main__":
    main(sys.argv)
