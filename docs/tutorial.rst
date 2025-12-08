Tutorial
========

This tutorial is intended to help you get started with the software, providing 
step-by-step instructions on how to input your data, run simulations, and evaluate your results.
If looking for specific topics, please refer to the relevant sections in the documentation.

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

To initiate a simulation using our software, you will need to provide an input file in JSON format. 
This file specifies all necessary parameters for the simulation, as well as the paths 
to additional input files such as imaging data or specific points to be analyzed. Below is a minimal example:

.. code-block:: json

  {
    "Electrodes": [
      {
        "Name": "BostonScientificVerciseDirected",
        "TipPosition": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0 },
        "Contacts": [
          {
            "Contact_ID": 1,
            "Active": true,
            "Voltage[V]": 1.0
          }
        ]
      }
    ],
    "MaterialDistribution": {
      "MRIPath": "path/to/mri.nii",
      "DiffusionTensorActive": true,
      "DTIPath": "path/to/dti.nii"
    },
    "OutputPath": "path/to/output",
    "StimulationSignal": {
      "Type": "Rectangle",
      "Frequency[Hz]": 130.0,
      "PulseWidth[us]": 60.0
    }
  }

This example demonstrates the minimal setup required to run a simulation. You can expand this JSON to include additional parameters as needed.

Provide Imaging Data
--------------------

Our software supports the use of segmented MRI and DTI images in Nifti format to account for the non-uniform distribution of brain tissue. 
These files can be stored anywhere, but their paths must be specified in the input JSON. If either file is not provided, the software will assume isotropic or homogeneous tissue properties, respectively.
Here is how to specify the paths in the JSON:

.. code-block:: bash

  "MaterialDistribution": {
    "MRIPath": "",
    "DiffusionTensorActive": false,
    "DTIPath": ""
  }

Starting simulation
--------------------

To start the simulation, run the software with the input file as an argument:

.. code-block:: bash

    $ ossdbs /absolute/path/to/input_file.json

The simulation will begin, and progress updates will be displayed in the console.
Once complete, the output files will be saved in the specified directory.

Evaluating Simulation Results
------------------------------

To review and analyze your simulation results, navigate to the output folder specified in the input JSON.
This folder contains all outputs, which can be imported into other analysis programs for further processing.
