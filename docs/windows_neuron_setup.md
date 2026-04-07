# NEURON Setup on Windows

This note explains how to make `NEURON` available to `OSS-DBS` on Windows for pathway activation modelling (PAM).
Following the steps below after you have successfully installed `OSS-DBS`.

## Download NEURON

Windows users can download the official installer from:

- <https://github.com/neuronsimulator/nrn/releases>

Use the Windows `.exe` installer from the release you want.

## Typical environment setup

After activating your Python environment, set the NEURON-related environment variables.

Example:

```powershell
$env:NEURONHOME = "C:\path\to\nrn"
$env:PYTHONPATH = "C:\path\to\nrn\lib\python"
$env:PATH = "C:\path\to\nrn\bin;$env:PATH"
```

Replace `C:\path\to\nrn` with the folder where NEURON was installed on your machine.

## Verify that NEURON works

Run:

```powershell
python -c "import neuron; print(neuron.__file__)"
```

If this succeeds, the Python bindings are available in the current environment.

You can also check whether `OSS-DBS` will detect PAM:

```powershell
python -c "import importlib.util; print(importlib.util.find_spec('neuron'))"
```

If the output is not `None`, `OSS-DBS` should see NEURON.