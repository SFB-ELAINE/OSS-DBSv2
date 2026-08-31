[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/SFB-ELAINE/OSS-DBSv2/main.svg)](https://results.pre-commit.ci/latest/github/SFB-ELAINE/OSS-DBSv2/main)

OSS-DBSv2
=========

This is the development version of the OSS-DBS rewrite.
Use with caution and please wait for the first official release before deploying it.
Bug reports are highly welcome, though!


Installation
------------

OSS-DBS is tested with Python 3.10, 3.11, 3.12, and 3.13.

**Windows users: Please install NEURON separately before
installing OSS-DBS. The instructions can be found [here](https://www.neuron.yale.edu/neuron/download).**


All other users can run

```
pip install ossdbs
```


## Developers

We recommend [`uv`](https://docs.astral.sh/uv/) for development. It is a fast,
drop-in replacement for `pip`/`virtualenv` that manages the virtual environment
and Python interpreter for you and uses the checked-in `uv.lock` for fully
reproducible installs. See the [installation instructions](https://docs.astral.sh/uv/getting-started/installation/).

Clone OSS-DBS into a local directory, `cd` into it, and create the environment
with the exact locked dependencies:

```
uv sync
```

`uv sync` creates a `.venv/` and installs OSS-DBS in editable mode. Select
optional dependency groups with `--extra`:

```
uv sync --extra test    # to run the test suite
uv sync --extra dev     # to develop OSS-DBS
uv sync --extra doc     # to locally build the docs
uv sync --extra all     # everything of the above
```

Run commands inside the environment with `uv run`, e.g. `uv run pytest` or
`uv run ossdbs input.json` (no manual activation needed).

<details>
<summary>Prefer plain <code>pip</code>? Click to expand.</summary>

Everything above also works with `pip`. Clone OSS-DBS, `cd` into the directory
and run one of:

```
pip install -e .            # base install
pip install -e ".[test]"    # to run the test suite
pip install -e ".[dev]"     # to develop OSS-DBS
pip install -e ".[doc]"     # to locally build the docs
pip install -e ".[all]"     # everything of the above
```

Note that `pip` resolves dependencies fresh and does not use `uv.lock`.
</details>

Run OSS-DBS
-----------

To run OSS-DBS, `cd` into the `input_files` directory, insert your parameters in the `input.json` 
and start the simulation with

```
ossdbs input.json
```

Also check out the `examples` directory and the documentation.

Development
-----------

The code development follows different coding styles that are checked
by git pre-commit hooks.
After `uv sync --extra dev` (or `pip install -e ".[dev]"`), activate the hooks
with `uv run pre-commit install` (or `pre-commit install`).

### Testing

OSS-DBS has two types of tests:

**Unit Tests** (`tests/`): Fast, isolated tests that run on every push/PR.
```bash
pytest                    # Run all unit tests
pytest --cov              # With coverage
```

**Simulation Tests** (`input_test_cases/`): Full simulation tests that validate the complete pipeline. These run only on pull requests to avoid slowing down regular development.
```bash
pytest input_test_cases/test_simulations.py -m simulation      # All simulation tests
pytest input_test_cases/test_simulations.py -m "simulation and not slow"  # Fast only
```

See [`input_test_cases/README.md`](input_test_cases/README.md) for detailed documentation on the simulation test suite, including how to add new tests.

### Directory Structure

- `tests/` - Unit tests (run on every commit)
- `input_test_cases/` - Simulation tests (run on PR only)
- `input_files/` - Reference input files
- `examples/` - Example code for users (may be resource-intensive)

Standalone executables
----------------------

`build_merged.py` freezes OSS-DBS and all of its dependencies into a
self-contained bundle with [PyInstaller](https://pyinstaller.org/), so that
end users do not need a Python toolchain. On Windows this is also how NEURON is
shipped: install and compile NEURON into the build environment first, then
freeze it into the bundle so users get PAM support without installing NEURON
themselves.

PyInstaller freezes whatever is present in the *current* environment, so
provision that environment reproducibly with `uv` and the checked-in
`uv.lock` before building:

```bash
uv sync                       # install the exact locked dependency graph
# (Windows only) install + compile NEURON into .venv here, see
# docs/windows_neuron_setup.md
uv run pip install pyinstaller
uv run python build_merged.py
```

The bundle is written to `dist/ossdbs_bundle/`. Because the environment comes
from `uv.lock`, the frozen dependency set is identical across rebuilds; run
`uv lock` (and commit the result) whenever you intentionally update a
dependency.
