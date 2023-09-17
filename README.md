OSS-DBSv2
=========

This is a first draft of the OSS-DBS rebuild.


Installation
------------

OSS-DBS is tested with Python 3.8, 3.9 and 3.10.

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
