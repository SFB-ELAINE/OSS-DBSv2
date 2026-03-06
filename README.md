[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/SFB-ELAINE/OSS-DBSv2/main.svg)](https://results.pre-commit.ci/latest/github/SFB-ELAINE/OSS-DBSv2/main)

OSS-DBSv2
=========

This is the development version of the OSS-DBS rewrite.
Use with caution and please wait for the first official release before deploying it.
Bug reports are highly welcome, though!


Installation
------------

OSS-DBS is tested with Python 3.8, 3.9, 3.10, and 3.11.

**Windows users and Mac users with Python 3.8: Please install NEURON separately before
installing OSS-DBS. The instructions can be found [here](https://www.neuron.yale.edu/neuron/download).**


All other users can run

```
pip install ossdbs
```


## Developers
To install OSS-DBS, clone it into a local directory,
`cd` into this directory and run

```
pip install -e .
```

To also run the test suite of OSS-DBS, run

```
pip install -e ".[test]"
```

To develop OSS-DBS, run

```
pip install -e ".[dev]"
```

To locally build the docs of OSS-DBS, run

```
pip install -e ".[doc]"
```

To do everything of the above, run

```
pip install -e ".[all]"
```

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
Install `pre-commit` via `pip install pre-commit` and run
`pre-commit install` to activate it.

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
