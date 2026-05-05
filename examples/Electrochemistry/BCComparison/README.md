# Boundary-condition comparison study

Six runs in two panels showing how electrode/contact boundary-condition
choices change the computed VTA. Geometry, mesh and signal are
identical across all six — only the contact BCs differ.

## Setup shared across all runs

- Electrode: `BostonScientificVerciseDirected`, geometry from
  `examples/ConvergenceStudy/VTA/`.
- Active stimulation: bipolar between segmented contacts —
  C2 at 0 V, C4 at +1 V (single 10 kHz frequency, voltage-controlled,
  QS mode).
- Mesh: HP refinement (Levels=2, Factor=0.125) on the default base
  mesh — convergence-study strategy XV (~60–120 k DOFs).
- Lattice for VTA: 60×60×240 points at 0.125 mm, threshold 200 V/m.

## Panels

| Panel | Run id          | Active C4 BC       | Non-active contacts (C1, C3, C5–C8) |
|-------|-----------------|--------------------|--------------------------------------|
| A     | `A1_floating`   | Dirichlet 1 V      | floating, no surface impedance       |
| A     | `A2_insulating` | Dirichlet 1 V      | insulating (Neumann)                 |
| A     | `A3_mild_Z`     | Dirichlet 1 V      | floating + R = 1 kΩ                  |
| B     | `B1_dirichlet`  | Dirichlet 1 V      | insulating                           |
| B     | `B2_mild_Z`     | 1 V + R = 1 kΩ     | insulating                           |
| B     | `B3_strong_Z`   | 1 V + R = 10 kΩ    | insulating                           |

C2 is the cathode (0 V) in every run.

`R` is the **specific** surface impedance (Ω·mm² in OSS-DBS's
convention; the code multiplies by contact area to get the Robin BC).

## Running

All six runs must execute sequentially (per repository CLAUDE.md
rule).

```bash
source venv/bin/activate
cd examples/Electrochemistry/BCComparison
bash run_all.sh
```

This runs Panel A, then Panel B, then aggregates results into
`bc_comparison_summary.csv`. Each run writes to `Results_<run_id>/`
with the usual OSS-DBS outputs (VCM report, lattice JSON, voxel
VTA NIfTI, impedance CSV, VTU export for ParaView).

## Quick-look results

`summarise.py` reads each `Results_*/` and writes
`bc_comparison_summary.csv` with columns `panel, run_id, dofs,
elements, time_s, ngs_vta_volume_mm3, 
impedance_real_Ohm, impedance_imag_Ohm`. It also prints the table to
stdout.

Sanity expectations before claiming success:

- Panel A: VTA volume ordering `A1_floating < A3_mild_Z < A2_insulating`
  (floating shunts current away from the active contact, insulating
  doesn't); impedance ordering reversed.
- Panel B: VTA volume ordering `B1_dirichlet > B2_mild_Z > B3_strong_Z`
  (more interface impedance → more voltage drops across the interface
  → smaller field in tissue); impedance ordering matches the same
  direction.

## Rendering the slice figures

Two parameterised renderers produce the same slice through the lead
(plane orthogonal to the lead axis, through the electrode centre)
with a shared camera, fixed colour bar and identical scalar ranges
across all six runs. Outputs go to `Figures_<run_id>/`.

| Script                       | Engine        | Output suffix       |
|------------------------------|---------------|---------------------|
| `render_figures.py`          | ParaView pvbatch | `*_slice.png`    |
| `render_figures_pyvista.py`  | PyVista (VTK) | `*_slice_pv.png`    |

A tissue-distribution slice (Gray matter / White matter / CSF colour
map of `material.vtu`) is rendered once into `Figures_material/` and
shared across all six runs since the geometry is identical.

### Slice + camera are derived from the JSON

`generate_paraview_camera.py` reads `base_settings.json` and writes
`paraview_camera.py` with slice plane and camera constants computed
from `Electrodes[0].Direction`, `TipPosition`, and
`BrainRegion.Center`. Both render scripts import that module, so
editing the JSON only requires a regenerate:

```bash
python generate_paraview_camera.py                       # default anchor = brain centre
python generate_paraview_camera.py --lead-offset-mm 3.25 # anchor on the contact array
python generate_paraview_camera.py --parallel-scale 18   # zoom out
```

OSS-DBS itself does not import ParaView; the generator only depends
on numpy + stdlib and emits Python text.

Both use:

- Slice origin = electrode centre `(11.6939, -15.9992, -10.1624)`,
  normal `(0, 0.805, -0.5101)` — orthogonal to the lead axis.
- Parallel projection, view-up = lead direction
  `(0.3027, 0.5101, 0.805)` so the lead is vertical on screen,
  `parallel_scale = 10 mm`.
- Potential range `[0, 1] V`, |E| range `[0, 300] V/m`.
- |E| converted from V/mm (FEM mm geometry) → V/m via ×1000.

### PyVista (in-venv, no ParaView needed)

```bash
source venv/bin/activate
cd examples/Electrochemistry/BCComparison
python render_figures_pyvista.py            # all six runs
python render_figures_pyvista.py A1_floating B1_dirichlet
```

### pvbatch (ParaView 6.1.0)

ParaView 5.12/5.13 OSMesa-MPI tarballs ship a broken
`libospray_module_ispc.so` that segfaults pvbatch on render. Use
**ParaView 6.1.0** (Qt6 build) and the offscreen flags below — this
is the only combination that actually works headless on this branch:

```bash
~/opt/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch \
    --force-offscreen-rendering --disable-xdisplay-test \
    render_figures.py
```

Add run ids as positional args to render a subset. The
`%-#6.3g` deprecation warnings on the scalar bar are cosmetic — they
will be silenced when ParaView upgrades the format string to
`std::format`.

## Final figure layout

Two rows × three columns; each cell pairs a potential slice and an
|E| slice for one run. Camera, colour preset and ranges are identical
across both panels so cells are directly comparable.
