NEURON install requires PATH modulation!

With the virtual environment activated, add NEURON's `bin` directory
(it contains `nrnivmodl` and `nrniv`) to your `PATH`:
```
export PATH="$(python -c 'import neuron, os; print(os.path.join(os.path.dirname(neuron.__file__), ".data", "bin"))')":$PATH
```
This works for any Python version. For example, for a virtualenv in your
home directory with Python 3.11 it expands to
`~/venv/lib/python3.11/site-packages/neuron/.data/bin`.
