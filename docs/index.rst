Welcome to OSS-DBSv2
====================

OSS-DBSv2 is an open-source simulation toolbox for deep brain stimulation (DBS).
It is designed for researchers and developers who want to model electric fields,
stimulation volumes, and pathway activation in patient-specific or experimental
setups.

The software can be used in two main ways:

- through `Lead-DBS <lead_dbs.html>`_ for a GUI-driven clinical and research workflow
- as a standalone Python and command-line toolbox for custom simulations,
  scripting, parameter sweeps, and method development

For most new users, the simplest path is to start with the standalone tutorial
or, if you already work in Lead-DBS, with the Lead-DBS integration page.

What OSS-DBSv2 can be used for
------------------------------

OSS-DBSv2 supports several common DBS modeling tasks:

- building model geometries from a simplified brain region or imaging data
- placing predefined or custom electrode models in the simulation domain
- assigning dielectric properties to tissue compartments
- solving the volume conductor problem in isotropic or anisotropic media
- evaluating fields on lattices, voxel grids, or pathway trajectories
- integrating with Lead-DBS for stimulation volume and pathway activation studies

Typical inputs
--------------

A simulation is typically defined by:

- a JSON input file describing geometry, electrodes, solver settings, and outputs
- segmented MRI data, optionally combined with DTI data
- a stimulation configuration with active contacts and signal settings
- an output directory for logs, meshes, and simulation results

Typical outputs
---------------

Depending on the workflow, OSS-DBSv2 can produce:

- electric field and potential solutions
- impedance estimates
- exported electrode and field data for visualization
- lattice- or pathway-based post-processing results
- log files and status files for reproducible batch runs

Where to start
--------------

If you are new to the project, the recommended path is:

1. :doc:`installation`
2. :doc:`tutorial`
3. :doc:`input_settings`
4. :doc:`examples`

If you want to drive OSS-DBSv2 from Python scripts instead of the CLI, see
:doc:`python_api`.

If you already use Lead-DBS, go directly to :doc:`lead_dbs`.

Documentation structure
-----------------------

The documentation is split into a user-oriented guide and a set of subsystem
reference pages. Start with the user guide if your goal is to run simulations.
Use the reference pages when you want more detail about geometry, materials,
signals, point analysis, or the solver stack.

For more details about the first version of OSS-DBS, see [Butenko2019]_.

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   installation
   input_settings
   lead_dbs
   tutorial
   examples
   python_api

.. toctree::
   :maxdepth: 1
   :caption: Documentation

   brain_geometry
   materials
   stimulation_signals
   electrodes
   point_analysis
   axon_models
   volume_conductor_model
   utils

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`

.. [Butenko2019] K. Butenko, C. Bahls, M. Schröder, R. Köhling and U. van Rienen,
   OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling,
   PLoS Comput Biol 16(7): e1008023.
   https://doi.org/10.1371/journal.pcbi.1008023
