"""SVG variant of the combined mesh-refinement overview with editable labels.

Same figure as ``render_mesh_refinement_combined.py`` (four refinement
panels + an hp magnifier inset), but the PyVista renders are kept as an
embedded raster while every text label is emitted as a real SVG ``<text>``
element via matplotlib (``svg.fonttype = 'none'``). That keeps the fonts
editable in Inkscape/Illustrator without re-rendering the meshes.

Reuses the geometry/camera constants and the ``_slab`` helper from
``render_mesh_refinement_combined`` so the two figures stay in sync.

Run generate_meshes.py, generate_fine_mesh.py and generate_hp_mesh.py first.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
sys.path.insert(0, _HERE)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pyvista as pv  # noqa: E402
import render_mesh_refinement as base  # noqa: E402
import render_mesh_refinement_combined as comb  # noqa: E402
from PIL import Image, ImageDraw  # noqa: E402

# keep text as editable <text> elements (system fonts), not outlined paths
plt.rcParams["svg.fonttype"] = "none"

OUTFILE = "mesh_refinement_overview_combined"


def _render_overview_notext(path):
    """Single-row whole-lead overview, no baked-in text labels."""
    lead = base.build_lead_parts()
    p = pv.Plotter(
        off_screen=True,
        shape=(1, len(comb.PANELS)),
        window_size=list(comb.WINDOW),
        border=False,
    )
    for idx, (result_dir, _label) in enumerate(comb.PANELS):
        p.subplot(0, idx)
        sl = comb._slab(result_dir)
        p.add_mesh(
            sl,
            scalars="material_real",
            cmap=base.CMAP,
            clim=base.CLIM,
            show_edges=True,
            edge_color=base.EDGE_COLOR,
            line_width=base.LINE_WIDTH,
            show_scalar_bar=False,
            **comb.LIGHT_KW,
        )
        base.add_lead(p, lead)
        p.set_background("white")
        p.enable_parallel_projection()
        p.camera_position = [
            (base.FOCAL[0], base.FOCAL[1] - 100.0, base.FOCAL[2]),
            base.FOCAL,
            (0.0, 0.0, 1.0),
        ]
        p.camera.parallel_scale = base.PARALLEL_SCALE
    p.screenshot(path)


def _render_inset_notext(path):
    """Zoomed contact-edge hp view with the edge marker, no baked-in text."""
    marker = pv.Line(
        (comb.CORNER[0], base.ELECTRODE_Y, comb.CONTACT_Z[0]),
        (comb.CORNER[0], base.ELECTRODE_Y, comb.CONTACT_Z[1]),
    )
    p = pv.Plotter(off_screen=True, window_size=list(comb.INSET_WINDOW), border=False)
    sl = comb._slab(comb.INSET_RESULT_DIR)
    p.add_mesh(
        sl,
        scalars="material_real",
        cmap=base.CMAP,
        clim=base.CLIM,
        show_edges=True,
        edge_color="#1a1a1a",
        line_width=1.4,
        show_scalar_bar=False,
        **comb.LIGHT_KW,
    )
    p.add_mesh(marker, color="#e8820c", line_width=7)
    p.set_background("white")
    p.enable_parallel_projection()
    p.camera_position = [
        (comb.CORNER[0], base.ELECTRODE_Y - 100.0, comb.CORNER[2]),
        comb.CORNER,
        (0.0, 0.0, 1.0),
    ]
    p.camera.parallel_scale = comb.INSET_SCALE
    p.screenshot(path)


def main():
    """Render text-free panels, composite the magnifier, add SVG text."""
    needed = [d for d, _ in comb.PANELS] + [comb.INSET_RESULT_DIR]
    missing = [d for d in needed if not os.path.exists(f"{d}/material.vtu")]
    if missing:
        raise SystemExit(f"Missing meshes (run the generate_*.py scripts): {missing}")

    _render_overview_notext("_overview_nt.png")
    _render_inset_notext("_inset_nt.png")

    img = Image.open("_overview_nt.png").convert("RGB")
    inset = Image.open("_inset_nt.png").convert("RGB")
    w, h = img.size
    pw = w / len(comb.PANELS)  # panel width in px
    half_h = base.PARALLEL_SCALE
    half_w = half_h * pw / h
    fx, fz = base.FOCAL[0], base.FOCAL[2]

    def to_px(x, z):
        px = comb.INSET_PANEL * pw + (x - (fx - half_w)) / (2 * half_w) * pw
        py = ((fz + half_h) - z) / (2 * half_h) * h
        return px, py

    bx0, by0 = to_px(comb.CORNER[0] - comb.INSET_HALF, comb.CORNER[2] + comb.INSET_HALF)
    bx1, by1 = to_px(comb.CORNER[0] + comb.INSET_HALF, comb.CORNER[2] - comb.INSET_HALF)

    draw = ImageDraw.Draw(img)
    draw.rectangle([bx0, by0, bx1, by1], outline=comb.ORANGE, width=3)  # source box

    iw, ih = inset.size
    pad = 14
    ix1 = int(comb.INSET_PANEL * pw + pw - pad)
    iy1 = h - pad
    ix0, iy0 = ix1 - iw, iy1 - ih
    draw.line([bx0, by1, ix0, iy0], fill=comb.ORANGE, width=2)
    draw.line([bx1, by1, ix1, iy0], fill=comb.ORANGE, width=2)
    img.paste(inset, (ix0, iy0))
    draw.rectangle([ix0, iy0, ix1 - 1, iy1 - 1], outline=comb.ORANGE, width=4)

    # raster mesh (with boxes/connectors) as the background; text as SVG <text>
    dpi = 100
    fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.imshow(img)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")

    for idx, (_rd, label) in enumerate(comb.PANELS):
        ax.text(
            idx * pw + 18,
            34,
            label,
            fontsize=20,
            color="black",
            va="top",
            ha="left",
            weight="bold",
        )
    ax.text(ix0 + 10, iy0 + 22, comb.INSET_LABEL, fontsize=14, color="black", va="top")
    ax.text(
        ix1 - 10,
        iy1 - 12,
        "contact edge",
        fontsize=13,
        color="#b5650a",
        va="bottom",
        ha="right",
    )

    fig.savefig(f"{OUTFILE}.svg")
    fig.savefig(f"{OUTFILE}_svg_preview.png", dpi=dpi)
    plt.close(fig)

    for tmp in ("_overview_nt.png", "_inset_nt.png"):
        os.remove(tmp)
    print(f"wrote {OUTFILE}.svg and {OUTFILE}_svg_preview.png")


if __name__ == "__main__":
    main()
