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

The code uses the `input_files` and `tests` directories to check the functionality
upon each commit. If you add a new feature, please add a test to `tests`.
Make sure that the test does not run long and does not consume much memory.
It shall be rather a unit test than a full simulation run.
Likewise, only change `input_files` after opening an issue.

The `examples` directory is meant for users to understand what has been implemented.
Place heavy and/or experimental code there (e.g., code that may crash or consume many resources).

The code development follows different coding styles that are checked
by git pre-commit hooks.
Install `pre-commit` via `pip install pre-commit` and run
`pre-commit install` to activate it. 
