# EasyVVUQ OSS-DBS Demo

This example shows a compact EasyVVUQ workflow around OSS-DBS impedance simulations. It varies three Cole–Cole parameters of the encapsulation layer:

- `sigma`: DC conductivity  
- `alpha_3`: third Cole–Cole dispersion exponent  
- `eps_delta_3`: third Cole–Cole dielectric increment  

The example uses `input_files/inputTest.json` as the base OSS-DBS input. The runner sets a thin encapsulation so the encapsulation-layer dielectric parameters affect the FEM problem.

---

## Files

- `EasyVVUQ_OSSDBS_Three_Parameter_Demo.ipynb`: PCE and SC EasyVVUQ workflow  
- `easyvvuq_input.template`: EasyVVUQ JSON input template  
- `uq_ossdbs_runner.py`: converts EasyVVUQ parameters to OSS-DBS input, runs OSS-DBS, and writes `output.csv`  

---

## Requirements

Use an environment where OSS-DBS runs, then install the temporary EasyVVUQ
compatibility branch:

```bash
python -m pip install "easyvvuq @ git+https://github.com/CheLamVien/EasyVVUQ.git@fix-pce-derivatives-numpoly2" pandas matplotlib scipy
```

This is the recommended installation for this example while EasyVVUQ's upstream
NumPy 2 support is still pending. The branch is related to the dependency update work in
[`UCL-CCS/EasyVVUQ#476`](https://github.com/UCL-CCS/EasyVVUQ/pull/476), but
keeps the changes scoped to dependency metadata and a PCE derivative fix. It
was tested with Python 3.13.13, NumPy 2.4.4, Chaospy 4.3.21, and Numpoly 1.3.9.

If you prefer the upstream EasyVVUQ 1.3 release, install it with NumPy 1.x:

```bash
python -m pip install "setuptools<81" "numpy<2" easyvvuq==1.3 chaospy pandas matplotlib scipy
```

EasyVVUQ 1.3 requires `numpy<2`, so this upstream-release path does not test the
new NumPy stack. The `setuptools<81` pin keeps `pkg_resources` available for the
Chaospy version resolved by EasyVVUQ 1.3.

---

## Known issues

### EasyVVUQ 1.3 bug

EasyVVUQ 1.3 contains a bug in:

```
easyvvuq/actions/execute_local.py
```

Replace:

```python
close(stdout)
close(stderr)
```

with:

```python
stdout.close()
stderr.close()
```

## Usage

Open and run:

```
EasyVVUQ_OSSDBS_Three_Parameter_Demo.ipynb
```
