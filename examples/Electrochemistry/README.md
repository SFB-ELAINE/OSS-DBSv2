# Electrochemistry: Surface Impedance at the Electrode-Tissue Interface

This example demonstrates the surface impedance (Robin boundary condition)
implementation in OSS-DBS, following the theory published in
[our previous paper](https://doi.org/10.1016/j.bioelechem.2023.108395).

## Background

In DBS, the electrode-tissue interface introduces an additional impedance
due to electrochemical processes (charge transfer, double-layer capacitance).
This is modeled as a surface impedance boundary condition (Robin BC) on the
active contact surfaces, rather than a simple Dirichlet BC.

The total macroscopic impedance measured between two contacts is:

```
Z_total = Z_tissue + Z_interface
```

where `Z_tissue` is the volume conductor impedance (from tissue conductivity)
and `Z_interface` is the electrode-tissue interface impedance (from the
electrochemical model, e.g., a constant-phase element).

## Input Files

### OSS-DBS JSON inputs

| File | Description |
|------|-------------|
| `input_no_interface.json` | No interface impedance (pure Dirichlet BC, contacts 1 & 8) |
| `input_no_interface_EQS.json` | No interface impedance (pure Dirichlet BC, contacts 1 & 8) |
| `input_interface_1kOhm.json` | 1 kOhm resistive interface impedance on contact 1 |
| `input_interface_1kOhm_EQS.json` | CPE interface (k=1e5, alpha=0.5) in EQS mode at 1 kHz |
| `input_interface_1kOhm_floating.json` | 1 kOhm interface on both floating contacts, current-controlled (2 mA) |
| `input_interface_1kOhm_EQS_floating.json` | CPE interface on floating contacts, EQS mode at 1 kHz |
| `input_interface_1kOhm_1active_1floating.json` | Mixed: 1 active contact + 1 floating with 1 kOhm interface |
| `input_interface_lempka2009.json` | Lempka 2009 CPE model (`CPE_dl`, k=1.5e6, alpha=0.8) at 100 Hz |
| `input_no_interface_lempka2009.json` | Same geometry, without surface impedance |
| `input_interface_lempka2009_spectrum.json` | Day 1 CPE (`dl_k=1.5e6, dl_alpha=0.8`) over 200 frequencies 10 Hz–100 kHz |
| `input_interface_lempka2009_day15_spectrum.json` | Day 15 CPE (`dl_k=1e6, dl_alpha=0.75`) over the same 200-frequency grid |
| `input_no_interface_lempka2009_spectrum.json` | Pure Dirichlet tissue-only reference over the same 200-frequency grid |
| `input_admittance_matrix.json` | 3-contact current-controlled setup for admittance/impedance matrix |

### Manual NGSolve scripts

These standalone scripts implement the same formulations directly in NGSolve,
serving as reference implementations for verification. Each script corresponds
to a JSON input and produces matching results (within ~1-2% due to mesh randomness).

| Script | JSON counterpart | Output CSV | Description |
|--------|-----------------|------------|-------------|
| `Simulation_no_Z.py` | `input_no_interface.json` | `results_no_Z_manual.csv` | Pure Dirichlet BC (no impedance) |
| `Simulation_hp.py` | `input_interface_1kOhm.json` | `results_manual.csv` | 1 kOhm Robin BC, voltage-controlled |
| `Simulation_hp_EQS.py` | `input_interface_1kOhm_EQS.json` | `results_Z_EQS.csv` | CPE Robin BC, EQS mode |
| `Simulation_hp_EQS_pure_float.py` | `input_interface_1kOhm_floating.json` | `results_Z_pure_float.csv` | Floating contacts + 1 kOhm Robin BC, current-controlled |

All scripts use HP-refinement (2 levels) on the active contact edges.

### Benchmark folder and regression check

The reference outputs of `Simulation_hp_EQS.py` and
`Simulation_hp_EQS_pure_float.py` are stored under `benchmark/`:

- `benchmark/results_Z_EQS.csv`
- `benchmark/results_Z_pure_float.csv`

To check the current implementation against these reference values, first
re-run the two manual simulations and the corresponding OSS-DBS JSON
inputs, then run `check_benchmarks.py`:

```bash
# Manual NGSolve scripts (produce CSVs in the current directory)
python Simulation_hp_EQS.py
python Simulation_hp_EQS_pure_float.py

# OSS-DBS JSON inputs (populate Results_*/ folders)
ossdbs input_interface_1kOhm_EQS.json
ossdbs input_interface_1kOhm_floating.json

# Compare everything against the benchmarks
python check_benchmarks.py
```

The check uses a 5% relative tolerance per physical quantity (`I`, `I1`,
`I2`, `V1`, `V2`, `Field`) to absorb mesh-randomness differences across
runs/machines. For the OSS-DBS comparisons, the check is sign-invariant
because OSS-DBS and the manual scripts use opposite sign conventions for
contact currents.

### Verification scripts

| Script | Purpose |
|--------|---------|
| `check_impedance.py` | Verifies the impedance decomposition `Z_total ≈ Z_tissue + Z_interface` for the OSS-DBS JSON examples (see *Verifying the impedance decomposition* below). |
| `check_benchmarks.py` | Compares both the manual NGSolve scripts and the OSS-DBS JSON runs against the reference CSVs in `benchmark/`. |
| `plot_impedance_spectrum.py` | Bode/Nyquist comparison of the full FEM run (`Results_LempkaImpedance_spectrum`) against `Z_tissue_FEM + Z_CPE_dl` for each Lempka 2009 parameter set (Day 1/15, pre/post stim). |
| `plot_impedance_spectrum_day1_vs_day15.py` | Overlays the Day 1 and Day 15 FEM spectra against their respective `Z_tissue_FEM + Z_CPE_dl` decompositions, plus a per-frequency relative-error plot. |

### Verification results

The following table compares OSS-DBS JSON results with the manual scripts
(different meshes, so results differ by ~1-2%):

| Scenario | OSS-DBS | Manual script |
|----------|---------|---------------|
| No interface | Z = 2174 Ohm | I = 0.466 mA (Z ~ 2146 Ohm) |
| 1 kOhm Robin | Z = 3216 Ohm | I = 0.314 mA (Z ~ 3185 Ohm) |
| CPE EQS | Z = 3110 - 896j Ohm | I = 0.299 + 0.087j mA (Z ~ 3083 - 897j Ohm) |
| Floating 1 kOhm | V_float = +/-4.26 V | V_float = +/-4.23 V |

## Verifying the impedance decomposition

To verify that `Z_total = Z_tissue + Z_interface`, run the Lempka example
with and without the surface impedance:

```bash
ossdbs input_interface_lempka2009.json
ossdbs input_no_interface_lempka2009.json
```

Then compare the results:

```python
import numpy as np
import impedancefitter as ifit

# Load results from impedance.csv files
Z_total  = ...  # from Results_LempkaImpedance/impedance.csv
Z_tissue = ...  # from Results_NoInterface_Lempka/impedance.csv

# Compute CPE impedance from impedancefitter
ecm = ifit.get_equivalent_circuit_model("CPE_dl")
Z_CPE = complex(ecm.eval(omega=2*np.pi*100, dl_k=1.5e6, dl_alpha=0.8))

# Check: Z_total ~ Z_tissue + Z_CPE (within ~1%)
Z_reconstructed = Z_tissue + Z_CPE
error = abs(Z_total - Z_reconstructed) / abs(Z_total)
```

The small (~1%) discrepancy arises because the Robin BC modifies the current
distribution in the tissue volume. The tissue impedance is not identical
between the two runs -- it shifts slightly when the surface impedance is
present, since the two components are coupled through the FEM solution.

For a frequency-swept version of the same check, run the spectrum inputs
and then one of the plot scripts:

```bash
ossdbs input_interface_lempka2009_spectrum.json
ossdbs input_interface_lempka2009_day15_spectrum.json
ossdbs input_no_interface_lempka2009_spectrum.json

python plot_impedance_spectrum.py             # Day 1 vs all 4 parameter sets
python plot_impedance_spectrum_day1_vs_day15.py  # Day 1 vs Day 15 overlay + error
```

Both plot scripts construct the decomposition as `Z_tissue_FEM(f) +
Z_CPE_dl(f)` using the FEM tissue spectrum (not an ECM tissue model) and
compare it against the full FEM spectrum. Max relative error is ~2% around
1–2 kHz for both Day 1 and Day 15.

## Admittance and impedance matrix for multicontact stimulation

For configurations with more than two contacts, the relationship between
voltages and currents is described by an admittance matrix Y (or its inverse,
the impedance matrix Z), rather than a single scalar impedance.

The file `input_admittance_matrix.json` demonstrates this with a 3-contact
current-controlled setup on a BostonScientificVercise electrode:

- E1C1: floating, +1 mA (source)
- E1C2: floating, -0.5 mA (sink)
- E1C7: active, 0V (ground), -0.5 mA (sink)

```bash
ossdbs input_admittance_matrix.json
```

This produces two CSV files in `Results_AdmittanceMatrix/`:

- `admittance_matrix.csv`: the full 3x3 admittance matrix Y (all contacts)
- `impedance_matrix.csv`: the reduced 2x2 impedance matrix Z (ground removed)

The admittance matrix Y is computed via the **superposition approach**: for N
contacts, N(N+1)/2 Dirichlet BVPs are solved with different voltage
configurations. Each entry Y_ij is derived from the power dissipated in the
corresponding BVP (see `examples/MulticontactCurrents/Superposition-approach.ipynb`
for the mathematical derivation).

Since the ground contact is included, the full Y matrix is singular (rows sum
to zero by current conservation). The impedance matrix is obtained by removing
the ground row/column from Y and inverting the resulting (N-1)x(N-1) matrix:

```
Z_reduced = inv(Y_reduced)
```

## How the surface impedance is implemented

The surface impedance is specified per contact in the JSON input:

```json
{
  "Contact_ID": 1,
  "Active": true,
  "Voltage[V]": 1.0,
  "SurfaceImpedance": {
    "Model": "CPE_dl",
    "Parameters": {"dl_k": 1.5e6, "dl_alpha": 0.8}
  }
}
```

The `Model` and `Parameters` fields are passed to
[impedancefitter](https://impedancefitter.readthedocs.io/) to evaluate the
specific impedance `Z_specific` (in Ohm) at each frequency. The surface
impedance used in the Robin BC is then `Z_surface = Z_specific * A_contact`
(in Ohm*mm^2), where `A_contact` is the contact area.

In the FEM formulation, the Robin BC adds the following terms:

- **Bilinear form:** `(1/Z_surface) * u * v * ds(contact)`
- **Linear form:** `(1/Z_surface) * V_contact * v * ds(contact)`

The total impedance computed by `compute_impedance()` includes both the
volume power dissipation and the interface power dissipation:

```
P_total = integral(conj(E) . J, volume) + integral((1/Z_s) |V_bc - u|^2, contact surface)
Z_total = |V_drop|^2 / P_total
```

## Reference

Zimmermann, J., Sahm, F., Arbeiter, N., Bathel, H., Song, Z., Bader, R.,
Jonitz-Heincke, A., van Rienen, U. (2023). Experimental and numerical
methods to ensure comprehensible and replicable alternating current electrical
stimulation experiments. *Bioelectrochemistry*, 151, 108395.
https://doi.org/10.1016/j.bioelechem.2023.108395
