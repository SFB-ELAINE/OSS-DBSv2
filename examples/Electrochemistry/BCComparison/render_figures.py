"""Render slice figures for the six BC-comparison runs with pvbatch.

Replaces the per-run Results_*.pvsm state files with a single
parameterised script. For every run id, a slice through the
lead is taken (plane orthogonal to the lead axis, through the
electrode centre) and two PNGs are produced:

    Figures_<run_id>/potential_slice.png
    Figures_<run_id>/E_field_slice.png

Camera, slice plane, colour map and scalar ranges are constants at
the top of this file so all six panels are directly comparable.

Run with:

    pvbatch render_figures.py                # all six runs
    pvbatch render_figures.py A1_floating    # one run
    pvbatch render_figures.py A1_floating B1_dirichlet
"""

import os
import sys

from paraview.simple import (
    Calculator,
    ColorBy,
    CreateRenderView,
    Delete,
    ExportView,
    GetColorTransferFunction,
    GetScalarBar,
    Hide,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    XMLUnstructuredGridReader,
)

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

# Fixed scalar ranges so all six panels share identical colour bars.
POTENTIAL_RANGE = (0.0, 1.0)  # Volts (Dirichlet BC = 1 V)
# E-field uses log scale, so the lower bound must be strictly positive.
E_FIELD_RANGE = (1.0, 300.0)  # V/m

# Material indices (see base_settings.json MRIMapping).
# 0 Unknown / outside, 1 Gray matter, 2 White matter, 3 CSF.
MATERIAL_RANGE = (-0.5, 3.5)
MATERIAL_RGB_POINTS = [
    -0.5,
    1.00,
    1.00,
    1.00,  # Unknown — white
    0.5,
    1.00,
    1.00,
    1.00,
    0.5,
    0.85,
    0.70,
    0.70,  # Gray matter — muted pink
    1.5,
    0.85,
    0.70,
    0.70,
    1.5,
    0.95,
    0.92,
    0.78,  # White matter — cream
    2.5,
    0.95,
    0.92,
    0.78,
    2.5,
    0.70,
    0.85,
    0.95,  # CSF — pale blue
    3.5,
    0.70,
    0.85,
    0.95,
]

# Output image size and presets.
IMAGE_SIZE = (1200, 1200)
# Diverging map for potential (signed BC), perceptually-uniform inferno
# for the strictly non-negative |E| magnitude.
POTENTIAL_PRESET = "Cool to Warm (Extended)"
E_FIELD_PRESET = "Inferno"


def _new_view():
    """Return a fresh render view configured to match the PVSM camera."""
    view = CreateRenderView()
    view.ViewSize = list(IMAGE_SIZE)
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    _apply_camera(view)
    return view


def _apply_camera(view):
    """(Re-)apply the shared camera. Called after Show() too because
    ParaView 6.x re-fits the camera to data extent on first display.
    """
    view.CameraParallelProjection = 1
    view.CameraPosition = list(CAMERA_POSITION)
    view.CameraFocalPoint = list(CAMERA_FOCAL_POINT)
    view.CameraViewUp = list(CAMERA_VIEW_UP)
    view.CameraParallelScale = CAMERA_PARALLEL_SCALE


def _slice(reader):
    """Build the cut filter at the fixed plane."""
    cut = Slice(Input=reader)
    cut.SliceType = "Plane"
    cut.SliceType.Origin = list(SLICE_ORIGIN)
    cut.SliceType.Normal = list(SLICE_NORMAL)
    cut.HyperTreeGridSlicer = "Plane"
    return cut


def _style_color_bar(view, field, preset, title, log_scale=False, fmt="{:.2f}"):
    """Apply consistent colour-bar styling, optionally on log scale."""
    lut = GetColorTransferFunction(field)
    # Reset to linear before applying the preset so that re-use of the global
    # LUT object across runs does not accumulate log-space transformations.
    lut.UseLogScale = 0
    lut.ApplyPreset(preset, True)
    if log_scale:
        lut.MapControlPointsToLogSpace()
        lut.UseLogScale = 1
    bar = GetScalarBar(lut, view)
    bar.Title = title
    bar.ComponentTitle = ""
    bar.LabelFormat = fmt
    bar.RangeLabelFormat = fmt
    bar.LabelColor = [1.0, 1.0, 1.0]
    bar.TitleColor = [1.0, 1.0, 1.0]
    bar.TitleFontSize = 32
    bar.LabelFontSize = 28
    bar.ScalarBarLength = 0.65
    bar.ScalarBarThickness = 30
    return lut, bar


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
    pot_reader = XMLUnstructuredGridReader(FileName=[pot_path])
    pot_reader.PointArrayStatus = ["potential_real"]
    pot_slice = _slice(pot_reader)

    view = _new_view()
    pot_disp = Show(pot_slice, view)
    pot_disp.Representation = "Surface"
    ColorBy(pot_disp, ("POINTS", "potential_real"))
    pot_disp.RescaleTransferFunctionToDataRange(False, True)
    lut, _ = _style_color_bar(
        view, "potential_real", POTENTIAL_PRESET, "Electric potential / V", fmt="{:.2f}"
    )
    lut.RescaleTransferFunction(*POTENTIAL_RANGE)
    pot_disp.SetScalarBarVisibility(view, True)
    Render(view)
    SaveScreenshot(
        os.path.join(out_dir, "potential_slice.png"),
        view,
        ImageResolution=list(IMAGE_SIZE),
        TransparentBackground=0,
    )
    ExportView(os.path.join(out_dir, "potential_slice.svg"), view=view)
    Hide(pot_slice, view)
    Delete(pot_slice)
    Delete(pot_reader)
    Delete(view)

    # ---- E-field magnitude slice (V/mm → V/m) ----
    e_reader = XMLUnstructuredGridReader(FileName=[e_path])
    e_reader.PointArrayStatus = ["E_field_real"]
    e_slice = _slice(e_reader)
    e_calc = Calculator(Input=e_slice)
    e_calc.AttributeType = "Point Data"
    e_calc.ResultArrayName = "E_Vm"
    e_calc.Function = "1000.0*mag(E_field_real)"

    view = _new_view()
    e_disp = Show(e_calc, view)
    e_disp.Representation = "Surface"
    ColorBy(e_disp, ("POINTS", "E_Vm"))
    e_disp.RescaleTransferFunctionToDataRange(False, True)
    lut, _ = _style_color_bar(
        view, "E_Vm", E_FIELD_PRESET, "E / V·m⁻¹", log_scale=True, fmt="{:.0f}"
    )
    lut.RescaleTransferFunction(*E_FIELD_RANGE)
    e_disp.SetScalarBarVisibility(view, True)
    Render(view)
    SaveScreenshot(
        os.path.join(out_dir, "E_field_slice.png"),
        view,
        ImageResolution=list(IMAGE_SIZE),
        TransparentBackground=0,
    )
    ExportView(os.path.join(out_dir, "E_field_slice.svg"), view=view)
    Hide(e_calc, view)
    Delete(e_calc)
    Delete(e_slice)
    Delete(e_reader)
    Delete(view)

    print(f"[ok]   {run_id}: wrote {out_dir}/{{potential,E_field}}_slice.{{png,svg}}")


def render_material(run_id: str = "A1_floating") -> None:
    """Render the tissue distribution slice once (identical across runs)."""
    mat_path = os.path.join(f"Results_{run_id}", "material.vtu")
    if not os.path.exists(mat_path):
        print(f"[skip] material: missing {mat_path}")
        return
    out_dir = "Figures_material"
    os.makedirs(out_dir, exist_ok=True)

    mat_reader = XMLUnstructuredGridReader(FileName=[mat_path])
    mat_reader.PointArrayStatus = ["material_real"]
    mat_slice = _slice(mat_reader)

    view = _new_view()
    mat_disp = Show(mat_slice, view)
    mat_disp.Representation = "Surface"
    ColorBy(mat_disp, ("POINTS", "material_real"))
    mat_disp.RescaleTransferFunctionToDataRange(False, True)
    lut = GetColorTransferFunction("material_real")
    lut.RGBPoints = list(MATERIAL_RGB_POINTS)
    lut.ColorSpace = "RGB"
    lut.RescaleTransferFunction(*MATERIAL_RANGE)
    bar = GetScalarBar(lut, view)
    bar.Title = "Material"
    bar.ComponentTitle = ""
    bar.LabelFormat = "{:<#6.0f}"
    bar.RangeLabelFormat = "{:<#6.0f}"
    bar.TitleFontSize = 32
    bar.LabelFontSize = 28
    bar.ScalarBarLength = 0.65
    bar.ScalarBarThickness = 30
    mat_disp.SetScalarBarVisibility(view, True)
    Render(view)
    SaveScreenshot(
        os.path.join(out_dir, "material_slice.png"),
        view,
        ImageResolution=list(IMAGE_SIZE),
        TransparentBackground=0,
    )
    ExportView(os.path.join(out_dir, "material_slice.svg"), view=view)
    Hide(mat_slice, view)
    Delete(mat_slice)
    Delete(mat_reader)
    Delete(view)

    print(f"[ok]   material: wrote {out_dir}/material_slice.{{png,svg}}")


def main(argv: list[str]) -> None:
    """Render every requested run; default to all six plus the material slice."""
    runs = argv[1:] if len(argv) > 1 else ALL_RUNS
    for run_id in runs:
        render_run(run_id)
    if len(argv) <= 1:
        render_material()


if __name__ == "__main__":
    main(sys.argv)
