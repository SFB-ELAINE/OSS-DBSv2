Getting Started
================

This tutorial is intended to help you get started with the software, providing 
step-by-step instructions on how to input your data, run simulations and evaluate your results.

Provide Input JSON 
-------------------

To initiate a simulation using our software, you will need to provide an input file in JSON format. 
This file will allow you to specify all necessary parameters for the simulation, as well as the path 
to additional input files such as imaging data or specific points to be analyzed. 
An standard input file without imaging data will have the following structure:

TODO: to be updated!

.. code-block:: bash

  "CaseGrounding": {
    "Active": false,
    "Current[A]": 0.0,
    "Voltage[V]": 0.0
  },
  "CurrentControlled": false,
  "DielectricModel": {"Type": ""
                    },
  "DiffusionTensorImage": {
    "Path": ""
  },
  "Electrodes": [
    {
      "Name": "BostonScientificVerciseDirected",
      "PathToCustomParameters": "",
      "Rotation[Degrees]": 0.0,
      "Direction": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0 },
      "TipPosition": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0 },

      "Contacts": [
        {
          "Contact_ID": 1,
          "Active": true,
          "Current[A]": 0.0,
          "Voltage[V]": 1.0,
          "Floating": false,
          "SurfaceImpedance[Ωm]": { "real": 0.0, "imag": 0.0 }
        }
      ]
    }
  ],
  "Contacts":{"MaxMeshSize": 0.1},
  "EncapsulationLayer": {
    "Thickness[mm]": 0.0,
    "Material": "",
    "MaxMeshSize": 0.0
  },
  "EQSMode": false,
  "Floating": {
    "Active": false,
    "FloatingImpedance": false
  },
  "MaterialDistribution": {
    "MRIPath": "",
    "DiffusionTensorActive": false,
    "DTIPath": ""
  },
  "Mesh": {
    "LoadMesh": false,
    "LoadPath": "",
    "MeshElementOrder": 0,
    "MeshingHypothesis": {
      "Type": "Default",
      "MaxMeshSize": 0.0
    },
    "SaveMesh": false
  },
  "OutputPath": "",
  "RegionOfInterest": {
    "Center": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0 },
    "Dimension": { "x[mm]": 10.0, "y[mm]": 10.0, "z[mm]": 10.0 }
  },
  "Solver": {
    "Type": "CG",
    "Preconditioner": "bddc",
    "PrintRates": false,
    "MaximumSteps": 10000,
    "Precision": 1e-12
  },
  "SpectrumMode": "FullSpectrum",
  "StimulationSignal": {
    "Type": "Rectangle",
    "Frequency[Hz]": 130.0,
    "PulseWidth[us]": 60.0,
    "PulseTopWidth[us]": 0.0,
    "CounterPulseWidth[us]": 0.0,
    "InterPulseWidth[us]": 0.0
  },
  "PointModel":{
    "Pathway": {
      "Active": false,
      "FileName": ""
    },
    "Lattice": {
      "Center": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0 },
      "Direction": { "x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0 },
      "PointDistance[mm]": 0.1,
      "Shape": { "x": 1, "y": 1, "z": 1 }
    }

Provide imaging data
--------------------

.. note::
    Providing DTI is not available yet.

In order to account for the non-uniform distribution of brain tissue, our software requires a segmented MRI of the brain in Nifti format. 
Additionally, a DTI image in Nifti format can be provided to incorporate the dispersive properties of the brain. 
These files can be stored anywhere, but the input JSON must include the file paths. 
If either one or both files are not provided, the software will assume isotropic or homogenous tissue properties respectivly.

.. code-block:: bash

  "MaterialDistribution": {
    "MRIPath": "",
    "DiffusionTensorActive": false,
    "DTIPath": ""
  }

Starting simulation
--------------------

To start the simulation, you can either navigate to the directory where the input JSON is located using the command 
line and run the software with the input file as an argument.

.. code-block:: bash

    $ ossdbs input_file.json 

Also it is possible to start the simulation form everywhere by giving the absolut path to the input JSON.

.. code-block:: bash

    $ ossdbs <path_to_input>/input_file.json 


The simulation will begin and progress updates will be displayed in the console. 
Once the simulation is complete, the output files will be saved in the specified directory in the input JSON. 


Evaluating simulation results
------------------------------

To review and analyze your simulation results, simply navigate to the output folder specified in the input JSON. 
This is where all outputs are stored, and the user can easily access and evaluate them. 
From this folder, you can import the data into other analysis programs for further processing.
