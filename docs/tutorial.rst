Tutorial
========

This tutorial is intended to help getting started with OSS-DBSv2.
It provides step-by-step instructions on how to set up your data, run simulations, and evaluate the results.
If you are looking for a specific topic, please refer to the relevant documentation sections below.

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

Getting Started
---------------

To run simulations using OSS-DBSv2 an input file in JSON format must be provided.
In this file, stimulation parameters and paths for input and output can be specified.
Using the example file ``inputTest.json`` from the ``input_files`` folder, a simulation can be run as follows:

.. code-block:: bash

    $ ossdbs /absolute/path/to/inputTest.json

Progress updates will be displayed in the console, and output and log files are stored in the specified directory.

OSS-DBSv2 can be used standalone or directly with reconstruction results from Lead-DBS.

Creating an Input File
----------------------

This section defines the brain region of interest for the simulation.

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

This section defines the electrode and contact properties.

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

Provide Imaging Data
--------------------

Segmented MRI and normalized DTI data in NIfTI format can be used to account for the non-uniform distribution of brain tissue.
These files can be stored anywhere, but their paths must be specified in the input JSON.
Additionally, the dielectric model can be specified to define the electrical properties of the tissue.

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

Stimulation Signal
------------------

The stimulation signal used in the simulation can be defined as follows:

.. code-block:: json

    "StimulationSignal": {
      "Type": "Rectangle",
      "Frequency[Hz]": 130.0,
      "PulseWidth[us]": 60.0,
      "PulseTopWidth[us]": 0.0,
      "CounterPulseWidth[us]": 0.0,
      "InterPulseWidth[us]": 0.0,
      "SpectrumMode": "OctaveBand",
      "CounterAmplitude": 1.0,
      "CutoffFrequency": 250000.0,
      "CurrentControlled": true
    },

Further Simulation Parameters
-----------------------------

More technical options, such as meshing and solver settings, can be adjusted as follows:

.. code-block:: json

    "Mesh": {
      "LoadMesh": false,
      "MeshElementOrder": 2,
      "MeshingHypothesis": {
        "Type": "Fine",
        "MaxMeshSize": 100
      },
      "SaveMesh": false
    },
    "EQSMode": false,
    "FEMOrder": 2,
    "Solver": {
      "Type": "CG",
      "Preconditioner": "bddc",
      "PreconditionerKwargs": {},
      "PrintRates": true,
      "MaximumSteps": 200,
      "Precision": 1e-9
    }

.. note::
   While the stimulation parameters and the paths to the imaging data are crucial, the technical parameters are safely covered by default settings.

Evaluating Simulation Results
------------------------------

To review and analyze your simulation results, navigate to the output folder specified in the input JSON.
Which outputs are generated can be configured as follows:

.. code-block:: json

    "OutputPath": "Results",
    "ComputeImpedance": true,
    "ExportVTK": true,
    "ExportElectrode": true

Additionally, the stimulation volume or pathway activation can be estimated using the implemented point models.
More details can be found in the :doc:`Point analysis <point_analysis>` and :doc:`Axon models <axon_models>` sections of the documentation.
