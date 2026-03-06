# OSS-DBS v2.0 Test Cases

This folder contains simulation test cases for validating OSS-DBS v2.0 functionality.

## Directory Structure

```
input_test_cases/
├── tests/                      # Test input configurations
│   ├── brain_material/         # Tissue property tests
│   ├── custom_parameters/      # Custom electrode/material tests
│   ├── case_grounding/         # Case grounding tests
│   ├── current_controlled/     # Current-controlled stimulation tests
│   ├── stimulation_signals/    # Signal processing tests
│   ├── floating_contacts/      # Floating conductor tests
│   ├── vta/                    # Volume of Tissue Activated tests
│   ├── pathway_activation/     # Pathway Activation Modeling tests
│   └── surface_impedance/      # Surface impedance tests
├── fixtures/                   # Shared test data
│   └── base_config.json        # Common default configuration
├── expected_outputs/           # Reference outputs for validation
├── manifest.yaml               # Test metadata and registry
├── conftest.py                 # Pytest configuration
├── test_simulations.py         # Pytest-based test runner
└── input_case*/                # Legacy test structure (deprecated)
```

## Running Tests

### Quick Start

```bash
# Run all simulation tests (on PR only by default)
pytest input_test_cases/test_simulations.py -m simulation

# Run only fast tests
pytest input_test_cases/test_simulations.py -m "simulation and not slow"

# Run specific category
pytest input_test_cases/test_simulations.py -k "brain_material"

# Run with verbose output
pytest input_test_cases/test_simulations.py -v --tb=short
```

### Test Markers

Tests are tagged with markers for selective execution:

| Marker | Description |
|--------|-------------|
| `simulation` | Full simulation tests (excluded from regular CI) |
| `slow` | Tests taking > 1 minute |
| `requires_neuron` | Tests requiring NEURON installation |
| `vta` | VTA computation tests |
| `pam` | Pathway Activation Modeling tests |
| `floating` | Floating contact tests |
| `surface_impedance` | Surface impedance tests |

### CI/CD Integration

- **Regular CI (push/PR)**: Runs unit tests only (`pytest` without simulation marker)
- **PR-only workflow**: Runs full simulation tests via `.github/workflows/simulation-tests.yml`

## Test Categories

### Brain Material (Case 1)
Tests homogeneous and inhomogeneous tissue models with ColeCole4 dielectric properties.
- Electrode: BostonScientificVercise
- Validates: Impedance computation

### Custom Parameters (Case 2)
Demonstrates custom electrode geometries and material models.
- Tests: Modified contact lengths, constant dielectric model
- Electrode: BostonScientificVerciseDirected with encapsulation

### Case Grounding (Case 3)
Tests monopolar stimulation with outer boundary grounding.
- Electrode: MicroProbesRodentElectrode
- Variants: QS mode, EQS mode

### Current Controlled (Case 4)
Tests current-controlled stimulation.
- Electrode: AbbottStJudeDirected6172
- Tests: Single contact, multi-contact current distribution

### Stimulation Signals (Case 5)
Tests frequency-domain signal processing.
- Signal: Rectangular pulse at 130 Hz
- Electrode: MedtronicSenSightB33005
- Uses OctaveBand spectrum approximation

### Floating Contacts (Case 6)
Tests floating conductor boundary conditions.
- Electrode: PINSMedicalL303
- Configuration: Active + floating + ground contacts

### VTA (Case 7)
Tests Volume of Tissue Activated computation.
- Electrode: DixiSEEG10
- Output: NIfTI format VTA
- Variants: Standard, out-of-core

### Pathway Activation (Case 8)
Tests neural pathway activation modeling.
- Electrode: Medtronic3387
- Requires: NEURON installation
- Input: HDF5 axon coordinates

### Surface Impedance (Case 9)
Tests FloatingImpedance formulation with electrode-tissue interface.
- Uses impedancefitter models (CPE, R)
- Mode: EQS with complex conductivity

## Adding New Tests

1. Create a JSON configuration in the appropriate `tests/<category>/` directory
2. Add test metadata to `manifest.yaml`
3. Add expected outputs to `expected_outputs/<category>/`
4. Run the test to verify: `pytest input_test_cases/test_simulations.py -k "your_test_id"`

### Minimal Test Configuration

Test files only need to specify what differs from `fixtures/base_config.json`:

```json
{
  "_comment": "Description of test",
  "_test_id": "unique_test_id",
  "Electrodes": [
    {
      "Name": "ElectrodeName",
      "Contacts": [
        {"Contact_ID": 1, "Active": true, "Voltage[V]": 1.0}
      ]
    }
  ]
}
```

## Migration from Legacy Structure

To migrate old `input_case*` directories:

```bash
python migrate_tests.py --dry-run  # Preview changes
python migrate_tests.py            # Execute migration
```

## Output Types

Simulation tests validate these output types:

| Output | File | Comparison Method |
|--------|------|-------------------|
| Impedance | `impedance.csv` | Numerical tolerance |
| VTA | `VTA_solution_Lattice.nii` | Dice coefficient |
| Floating Potentials | `floating_potentials.csv` | Numerical tolerance |

## Environment Variables

| Variable | Description |
|----------|-------------|
| `OSSDBS_KEEP_TEST_OUTPUTS` | Set to `true` to preserve test outputs |
