"""Combined mesh-refinement overview with an hp magnifier inset.

Single-row figure of four refinement strategies in canonical order
(Default I, Global fine II, Edge refinement V, Material refinement VIII).
hp refinement is not a full panel: it is only visible right at the contact
edge, so it is shown as a magnifier inset on the edge-refinement panel
(both refine the contact region), connected to its source box. The PyVista
renders are composited with PIL.

Run generate_meshes.py, generate_fine_mesh.py and generate_hp_mesh.py first.
"""

import os

import pyvista as pv
import render_mesh_refinement as base
from PIL import Image, ImageDraw

# (result dir, label) — each panel shows one representative mesh and a
# label spanning the canonical numerals of all strategies of that type:
#   I        default (baseline)
#   II-III   global fine / very fine
#   IV-VII   edge refinements (perimeter/20, /50, /75, voxel)
#   VIII-XII material refinements (material, edge+material combos)
#   XIII-XIV hp refinements (shown in the inset)
PANELS = [
    ("Results_PAM_default", "Default (I)"),
    ("Results_PAM_fine", "Global fine (II-III)"),
    ("Results_PAM_fine_edge_refinement", "Edge refinement (IV-VII)"),
    ("Results_PAM_material_refinement", "Material refinement (VIII-XII)"),
]
# hp is not a full panel; it is shown only as a magnifier inset on the
# edge-refinement panel (edge refinement and hp both refine the contacts)
INSET_PANEL = 2
INSET_RESULT_DIR = "Results_PAM_hp_refinement"
INSET_LABEL = "hp refinement (XIII-XIV)"
WINDOW = (3200, 950)

LEAD_R = 0.635
CORNER = (base.ELECTRODE_X + LEAD_R, base.ELECTRODE_Y, -7.02)  # a contact-rim corner
CONTACT_Z = (-8.52, -7.02)
INSET_WINDOW = (380, 400)
INSET_SCALE = 0.42
INSET_HALF = INSET_SCALE  # source-box half-size (mm) = inset half-height

ORANGE = (232, 130, 12)
OUTFILE = "mesh_refinement_overview_combined"

# Thin extracted-cell slab + soft lighting => 3-D shaded cells (ParaView's
# "extract cells by region"). SLAB_HALF small = one cell layer, little
# overlap; high ambient / low diffuse keeps the material colours readable.
SLAB_HALF = 0.1  # mm half-thickness
LIGHT_KW = {"lighting": True, "ambient": 0.6, "diffuse": 0.4, "specular": 0.0}


def _slab(result_dir):
    """Extract a thin slab of whole cells around the electrode plane.

    Two crinkle clips keep the tetrahedra whole instead of cutting them, so
    they render as shaded 3-D facets rather than a flat cross-section.
    """
    mesh = pv.read(os.path.join(result_dir, "material.vtu"))
    n = base.SLICE_NORMAL
    ox, oy, oz = base.SLICE_ORIGIN
    return mesh.clip(
        normal=n, origin=(ox, oy - SLAB_HALF, oz), invert=False, crinkle=True
    ).clip(normal=n, origin=(ox, oy + SLAB_HALF, oz), invert=True, crinkle=True)


def _render_overview(path):
    """Single-row whole-lead overview of all four strategies."""
    lead = base.build_lead_parts()
    p = pv.Plotter(
        off_screen=True, shape=(1, len(PANELS)), window_size=list(WINDOW), border=False
    )
    for idx, (result_dir, label) in enumerate(PANELS):
        p.subplot(0, idx)
        sl = _slab(result_dir)
        p.add_mesh(
            sl,
            scalars="material_real",
            cmap=base.CMAP,
            clim=base.CLIM,
            show_edges=True,
            edge_color=base.EDGE_COLOR,
            line_width=base.LINE_WIDTH,
            show_scalar_bar=False,
            **LIGHT_KW,
        )
        base.add_lead(p, lead)
        p.add_text(
            label, font_size=18, color="black", position="upper_left", shadow=True
        )
        p.set_background("white")
        p.enable_parallel_projection()
        p.camera_position = [
            (base.FOCAL[0], base.FOCAL[1] - 100.0, base.FOCAL[2]),
            base.FOCAL,
            (0.0, 0.0, 1.0),
        ]
        p.camera.parallel_scale = base.PARALLEL_SCALE
    p.screenshot(path)


def _render_inset(path):
    """Zoomed contact-edge view of the hp mesh (no lead), with an edge marker."""
    marker = pv.Line(
        (CORNER[0], base.ELECTRODE_Y, CONTACT_Z[0]),
        (CORNER[0], base.ELECTRODE_Y, CONTACT_Z[1]),
    )
    p = pv.Plotter(off_screen=True, window_size=list(INSET_WINDOW), border=False)
    sl = _slab(INSET_RESULT_DIR)
    p.add_mesh(
        sl,
        scalars="material_real",
        cmap=base.CMAP,
        clim=base.CLIM,
        show_edges=True,
        edge_color="#1a1a1a",
        line_width=1.4,
        show_scalar_bar=False,
        **LIGHT_KW,
    )
    p.add_mesh(marker, color="#e8820c", line_width=7)
    p.add_text(INSET_LABEL, font_size=13, color="black", position="upper_left")
    p.add_text("contact edge", font_size=12, color="#b5650a", position="lower_right")
    p.set_background("white")
    p.enable_parallel_projection()
    p.camera_position = [
        (CORNER[0], base.ELECTRODE_Y - 100.0, CORNER[2]),
        CORNER,
        (0.0, 0.0, 1.0),
    ]
    p.camera.parallel_scale = INSET_SCALE
    p.screenshot(path)


def main():
    """Render the overview + inset and composite the magnifier figure."""
    needed = [d for d, _ in PANELS] + [INSET_RESULT_DIR]
    missing = [d for d in needed if not os.path.exists(f"{d}/material.vtu")]
    if missing:
        raise SystemExit(f"Missing meshes (run the generate_*.py scripts): {missing}")

    _render_overview("_overview.png")
    _render_inset("_inset_raw.png")

    img = Image.open("_overview.png").convert("RGB")
    inset = Image.open("_inset_raw.png").convert("RGB")
    w, h = img.size
    pw = w / len(PANELS)  # panel width in px
    half_h = base.PARALLEL_SCALE
    half_w = half_h * pw / h
    fx, fz = base.FOCAL[0], base.FOCAL[2]

    def to_px(x, z):
        px = INSET_PANEL * pw + (x - (fx - half_w)) / (2 * half_w) * pw
        py = ((fz + half_h) - z) / (2 * half_h) * h
        return px, py

    bx0, by0 = to_px(CORNER[0] - INSET_HALF, CORNER[2] + INSET_HALF)
    bx1, by1 = to_px(CORNER[0] + INSET_HALF, CORNER[2] - INSET_HALF)

    draw = ImageDraw.Draw(img)
    draw.rectangle([bx0, by0, bx1, by1], outline=ORANGE, width=3)  # source box

    iw, ih = inset.size
    pad = 14
    ix1 = int(INSET_PANEL * pw + pw - pad)
    iy1 = h - pad
    ix0, iy0 = ix1 - iw, iy1 - ih
    # connectors: source-box bottom corners -> inset top corners
    draw.line([bx0, by1, ix0, iy0], fill=ORANGE, width=2)
    draw.line([bx1, by1, ix1, iy0], fill=ORANGE, width=2)
    img.paste(inset, (ix0, iy0))
    draw.rectangle([ix0, iy0, ix1 - 1, iy1 - 1], outline=ORANGE, width=4)

    img.save(f"{OUTFILE}.png")
    for tmp in ("_overview.png", "_inset_raw.png"):
        os.remove(tmp)
    print(f"wrote {OUTFILE}.png")


if __name__ == "__main__":
    main()
