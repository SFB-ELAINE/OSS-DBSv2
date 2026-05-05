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

Use an environment where OSS-DBS runs, then install dependencies:

```bash
python -m pip install easyvvuq==1.3 chaospy pandas matplotlib scipy
```

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