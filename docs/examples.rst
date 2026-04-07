.. _examples:

Examples
========

Different examples are provided with the software to help users move from a
first run to more specialized workflows. Some examples are simple scripts,
while others are notebooks or full JSON-based simulation setups.

Recommended starting points
---------------------------

If you are new to OSS-DBSv2, start here:

- :doc:`tutorial` for a first standalone simulation based on the shipped example
  input file
- :doc:`examples/example_JD` for a documented human DBS example
- :doc:`examples/example_rat` for a rodent example

Additional repository examples
------------------------------

The repository also contains example collections for specific subsystems:

- ``examples/BrainGeometryAPI`` for geometry construction
- ``examples/ElectrodesAPI`` for electrode placement and customization
- ``examples/DielectricModelAPI`` for tissue property models
- ``examples/ImageAPI`` for MRI and DTI handling
- ``examples/MeshAPI`` for mesh generation and refinement
- ``examples/StimulationSignalAPI`` for signal generation
- ``examples/MulticontactCurrents`` for multicontact stimulation studies
- ``examples/OptimizeSettings`` for exploratory optimization workflows

These examples are especially useful once the basic CLI workflow from the
tutorial is already familiar.

.. toctree::
   :maxdepth: 1

   examples/example_JD
   examples/example_rat
   examples/example_test_cases

.. note::
   More examples are stored in the ``examples`` and the ``input_test_cases`` directories within the repository.
