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

The ``examples/`` directory in the repository contains script and notebook
collections for specific subsystems. Each subdirectory is self-contained
and can be run after installation. The links below point to the
corresponding documentation pages for background.

- ``examples/BrainGeometryAPI`` — building brain regions with electrodes
  and encapsulation layers. See :doc:`brain_geometry`.
- ``examples/ElectrodesAPI`` — placing electrodes, including directed and
  custom models. See :doc:`electrodes`.
- ``examples/DielectricModelAPI`` — Cole-Cole and constant dielectric
  models, custom parameter overrides. See :doc:`materials`.
- ``examples/ImageAPI`` — MRI loading, voxel mapping, and DTI anisotropy
  (including the ``example_dti_mask.py`` masking demo). See :doc:`materials`
  (:ref:`dti-anisotropy`).
- ``examples/MeshAPI`` — mesh generation, local mesh sizes, HP refinement,
  and adaptive refinement. See :ref:`mesh-refinement`.
- ``examples/StimulationSignalAPI`` — rectangular and biphasic pulse
  generation, octave-band mode. See :doc:`stimulation_signals`.
- ``examples/MulticontactCurrents`` — current-controlled multicontact
  setups with floating and superposition approaches. See
  :ref:`stimulation-modes`.
- ``examples/OptimizeSettings`` — exploratory optimisation of stimulation
  parameters. See :doc:`python_api`.

.. toctree::
   :maxdepth: 1

   examples/example_JD
   examples/example_rat
   examples/example_test_cases

.. note::
   More examples are stored in the ``examples`` and the ``input_test_cases`` directories within the repository.
