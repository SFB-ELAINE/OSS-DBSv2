# Boundary-condition comparison study

Six runs in two panels showing how electrode/contact boundary-condition
choices change the computed VTA. Geometry, mesh and signal are
identical across all six — only the contact BCs differ.

See `PLAN.md` for the rationale and design decisions.

## Setup shared across all runs

- Electrode: `BostonScientificVerciseDirected`, geometry from
  `examples/ConvergenceStudy/VTA/`.
- Active stimulation: contact 2 at +1 V, brain-surface case grounded
  at 0 V (single 10 kHz frequency, voltage-controlled, QS mode).
- Mesh: HP refinement (Levels=2, Factor=0.125) on the default base
  mesh — convergence-study strategy XV (~60–120 k DOFs).
- Lattice for VTA: 60×60×240 points at 0.125 mm, threshold 200 V/m.

## Panels

| Panel | Run id          | Active C2 BC       | Non-active contacts (C1, C3–C8) |
|-------|-----------------|--------------------|----------------------------------|
| A     | `A1_floating`   | Dirichlet 1 V      | floating, no surface impedance   |
| A     | `A2_insulating` | Dirichlet 1 V      | insulating (Neumann)             |
| A     | `A3_mild_Z`     | Dirichlet 1 V      | floating + R = 1 kΩ              |
| B     | `B1_dirichlet`  | Dirichlet 1 V      | insulating                       |
| B     | `B2_mild_Z`     | 1 V + R = 1 kΩ     | insulating                       |
| B     | `B3_strong_Z`   | 1 V + R = 10 kΩ    | insulating                       |

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
elements, time_s, ngs_vta_volume_mm3, voxel_vta_volume_mm3,
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

## Final figure (ParaView)

Two rows × three columns; each cell shows the ±200 V/m VTA isosurface
plus the lead geometry from the VTU export, with consistent camera and
colour scale across both panels.
