Tutorial
========

This tutorial is intended to get a new user from installation to a first
successful simulation as quickly as possible. It focuses on the standalone
command-line workflow and uses the example files shipped with the repository.

Quickstart
----------

OSS-DBSv2 runs from a JSON input file. The repository already includes a
working example in ``input_files/inputTest.json``.

After installation, move into the example directory and start the simulation:

.. code-block:: bash

    $ cd input_files
    $ ossdbs inputTest.json

Progress information is printed to the console. Output files are written to the
directory specified by ``OutputPath`` in the JSON file. In the provided example,
that directory is ``Results`` inside ``input_files``.

If the run finishes successfully, you should see a results directory with log
files and exported simulation outputs. That is enough for a first sanity check
before you start changing inputs.

What the example run contains
-----------------------------

The example input file demonstrates a basic voltage-controlled simulation with:

- a spherical brain region centered around the electrode
- one Boston Scientific Vercise electrode
- MRI-based tissue labeling from ``Butenko_segmask.nii.gz``
- a Cole-Cole dielectric model
- a finite-element solve with configurable mesh and solver settings

This makes it a good starting point for understanding the structure of an
OSS-DBSv2 simulation before adapting it to new subjects or experiments.

Core input sections
-------------------

The JSON file is organized into a few main blocks. These are the most important
ones to inspect first.

Brain region
^^^^^^^^^^^^

This section defines the size and location of the simulated domain:

.. code-block:: json

    "BrainRegion": {
      "Center": {
        "x[mm]": -9.48,
        "y[mm]": 11.61,
        "z[mm]": 4.68
      },
      "Dimension": {
        "x[mm]": 40.0,
        "y[mm]": 40.0,
        "z[mm]": 40.0
      },
      "Shape": "Sphere"
    },

Electrode and contacts
^^^^^^^^^^^^^^^^^^^^^^

This section defines the electrode model, implantation direction, tip position,
and contact-level activation settings:

.. code-block:: json

    "Electrodes": [
      {
        "Name": "BostonScientificVercise",
        "Rotation[Degrees]": 0.0,
        "Direction": {
          "x[mm]": 0.0,
          "y[mm]": 0.0,
          "z[mm]": 1.0
        },
        "TipPosition": {
          "x[mm]": -9.48,
          "y[mm]": 11.61,
          "z[mm]": 4.68
        },
        "Contacts": [
          {
            "Contact_ID": 1,
            "Active": true,
            "Current[A]": 0.0,
            "Voltage[V]": 1.0,
            "Floating": false
          },
          {
            "Contact_ID": 2,
            "Active": true,
            "Current[A]": 0.0,
            "Voltage[V]": 0.0,
            "Floating": false
          }
        ]
      }
    ],

Imaging data and tissue mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Segmented MRI and optional DTI data in NIfTI format define the spatial material
distribution. The dielectric model specifies how tissue properties are computed.

.. code-block:: json

    "MaterialDistribution": {
      "MRIPath": "path/to/segmask.nii.gz",
      "MRIMapping": {
        "Unknown": 0,
        "CSF": 1,
        "White matter": 2,
        "Gray matter": 3,
        "Blood": 4
      },
      "DiffusionTensorActive": true,
      "DTIPath": "path/to/dti.nii.gz"
    },
    "DielectricModel": {
      "Type": "ColeCole4"
    },

Mesh, solver, and outputs
^^^^^^^^^^^^^^^^^^^^^^^^^

Numerical settings define how the problem is discretized and solved. The output
section determines where results are written and which exports are enabled.

Typical fields to look for are:

- ``Mesh`` for mesh generation and refinement
- ``Solver`` for the linear solver and preconditioner
- ``StimulationSignal`` for the excitation model
- ``OutputPath`` for results
- ``ExportVTK`` and ``ExportElectrode`` for additional exports

A stimulation signal block can look like this:

.. code-block:: json

    "StimulationSignal": {
      "Type": "Rectangle",
      "Frequency[Hz]": 130.0,
      "PulseWidth[us]": 60.0,
      "CurrentControlled": true
    },

At this stage, it is usually best to keep the default mesh and solver settings
and first confirm that your geometry, imaging paths, and stimulation settings
are correct.

Adapting the example to your own data
-------------------------------------

When moving from the example dataset to your own setup, the usual order is:

1. Replace ``MRIPath`` and, if needed, ``DTIPath``.
2. Update ``MRIMapping`` so label values match your segmentation.
3. Change the ``BrainRegion`` center and dimensions.
4. Replace the example electrode with the implanted model and correct position.
5. Adjust contact activation and stimulation settings.
6. Run the simulation once with conservative defaults before tuning the mesh or solver.

Standalone and Lead-DBS workflows
---------------------------------

OSS-DBSv2 can be used directly or through Lead-DBS:

- Use standalone mode when you want explicit control over settings, automation,
  scripting, or method development.
- Use Lead-DBS when you want image preprocessing, electrode reconstruction, and
  tight GUI integration.

If you already work in Lead-DBS, continue with :doc:`lead_dbs`.

Next steps
----------

After your first successful run, the most useful follow-up pages are:

- :doc:`input_settings`
- :doc:`examples`
- :doc:`electrodes`
- :doc:`volume_conductor_model`
- :doc:`python_api` — for scripting, parameter sweeps, and custom pipelines

Related Documentation
---------------------

- :doc:`Brain geometry <brain_geometry>`
- :doc:`Materials <materials>`
- :doc:`Stimulation signals <stimulation_signals>`
- :doc:`Electrodes <electrodes>`
- :doc:`Point analysis <point_analysis>`
- :doc:`Axon models <axon_models>`
- :doc:`Volume conductor model <volume_conductor_model>`
- :doc:`Utilities <utils>`
