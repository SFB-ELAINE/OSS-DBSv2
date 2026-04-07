# PAM Examples

This folder contains a self-contained example of pathway activation modelling that can be run directly from `examples/PAM/`.

The folder now includes the required input data:

- `Allocated_axons.h5`
- `Allocated_axons_parameters.json`
- `segmask.nii.gz`
- `oss_dbs_parameters.json`

These files were copied from the convergence-study PAM example. 


## Prerequisites

- a working OSS-DBS Python environment
- NEURON available in the active environment

For Windows-specific NEURON setup, see:

- `docs/windows_neuron_setup.md`

## Step 1: Run OSS-DBS on the pathway example

From the repository root:

```
cd examples\PAM
python -m ossdbs.main .\oss_dbs_parameters.json
```

This creates the pathway-sampled time-domain solution in:

- `examples/PAM/Results_rh/oss_time_result_PAM.h5`

## Step 2: Run the bundled PAM test

Use the generic test:

```
python .\NEURON_Test\test.py
```

Or use the model-specific McNeal test:

```
python .\NEURON_McNeal_Test\test.py
```

