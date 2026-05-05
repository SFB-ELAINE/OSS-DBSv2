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
"""

import os
import sys

import numpy as np
import pyvista as pv
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

POTENTIAL_RANGE = (0.0, 1.0)
# E-field uses log scale, so the lower bound must be strictly positive.
E_FIELD_RANGE = (1.0, 300.0)

# Material indices (see base_settings.json MRIMapping).
# 0 Unknown / outside, 1 Gray matter, 2 White matter, 3 CSF, 4 Blood.
# Blood is unused here (encapsulation thickness = 0).
MATERIAL_LABELS = ["Unknown", "Gray matter", "White matter", "CSF"]
MATERIAL_COLORS = [
    (1.00, 1.00, 1.00),  # Unknown — white
    (0.85, 0.70, 0.70),  # Gray matter — muted pink
    (0.95, 0.92, 0.78),  # White matter — cream
    (0.70, 0.85, 0.95),  # CSF — pale blue
]
MATERIAL_RANGE = (-0.5, 3.5)

IMAGE_SIZE = (1200, 1200)
# Diverging map for the (signed) potential — matplotlib's "coolwarm"
# matches ParaView's "Cool to Warm (Extended)" closely. Inferno for
# the strictly non-negative |E| magnitude (perceptually uniform,
# colour-blind friendly, identical preset name in ParaView).
POTENTIAL_CMAP = "coolwarm"
E_FIELD_CMAP = "inferno"

pv.OFF_SCREEN = True


def _setup_camera(plotter: pv.Plotter) -> None:
    """Apply the shared camera so all figures match."""
    plotter.enable_parallel_projection()
    plotter.camera_position = [
        list(CAMERA_POSITION),
        list(CAMERA_FOCAL_POINT),
        list(CAMERA_VIEW_UP),
    ]
    plotter.camera.parallel_scale = CAMERA_PARALLEL_SCALE
    plotter.background_color = "white"


def _slice(path: str) -> pv.PolyData:
    """Read a VTU and return the cut at the fixed plane."""
    mesh = pv.read(path)
    return mesh.slice(normal=SLICE_NORMAL, origin=SLICE_ORIGIN)


def _save(plotter: pv.Plotter, path: str) -> None:
    """Render and write a PNG, then close."""
    plotter.show(screenshot=path, auto_close=True)


def render_run(run_id: str) -> None:
    """Render potential and E-field slices for one run."""
    result_dir = f"Results_{run_id}"
    out_dir = f"Figures_{run_id}"
    os.makedirs(out_dir, exist_ok=True)

    pot_path = os.path.join(result_dir, "potential.vtu")
    e_path = os.path.join(result_dir, "E-field.vtu")
    if not (os.path.exists(pot_path) and os.path.exists(e_path)):
        print(f"[skip] {run_id}: missing potential.vtu or E-field.vtu")
        return

    # ---- potential slice ----
    pot_slice = _slice(pot_path)
    p = pv.Plotter(off_screen=True, window_size=list(IMAGE_SIZE))
    p.add_mesh(
        pot_slice,
        scalars="potential_real",
        cmap=POTENTIAL_CMAP,
        clim=POTENTIAL_RANGE,
        scalar_bar_args={
            "title": "Electric potential / V",
            "fmt": "%.2f",
            "n_labels": 5,
            "color": "white",
        },
    )
    _setup_camera(p)
    _save(p, os.path.join(out_dir, "potential_slice_pv.png"))

    # ---- E-field magnitude slice ----
    # FEM geometry is in mm so E_field_real is in V/mm; convert to V/m.
    e_slice = _slice(e_path)
    e_slice["|E|"] = 1000.0 * np.linalg.norm(e_slice["E_field_real"], axis=1)
    p = pv.Plotter(off_screen=True, window_size=list(IMAGE_SIZE))
    p.add_mesh(
        e_slice,
        scalars="|E|",
        cmap=E_FIELD_CMAP,
        clim=E_FIELD_RANGE,
        log_scale=True,
        scalar_bar_args={
            "title": "log(E) / V·m⁻¹",
            "fmt": "%.0f",
            "n_labels": 5,
            "color": "white",
        },
    )
    _setup_camera(p)
    _save(p, os.path.join(out_dir, "E_field_slice_pv.png"))

    print(f"[ok]   {run_id}: wrote {out_dir}/*_slice_pv.png")


def render_material(run_id: str = "A1_floating") -> None:
    """Render the tissue distribution slice once (identical across runs)."""
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
            "title": "Material",
            "n_labels": len(MATERIAL_LABELS),
            "fmt": "%.0f",
        },
    )
    _setup_camera(p)
    _save(p, os.path.join(out_dir, "material_slice_pv.png"))
    print(f"[ok]   material: wrote {out_dir}/material_slice_pv.png")


def main(argv: list[str]) -> None:
    """Render every requested run; default to all six plus the material slice."""
    runs = argv[1:] if len(argv) > 1 else ALL_RUNS
    for run_id in runs:
        render_run(run_id)
    if len(argv) <= 1:
        render_material()


if __name__ == "__main__":
    main(sys.argv)
