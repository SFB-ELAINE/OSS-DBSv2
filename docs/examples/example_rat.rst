Rodent studies
==============

In this example we show how to run a study in a rodent brain with OSS-DBS.
Here, we simulate the treatment of a wistar rat using a MRI based on Johnson et.al rat atlas.
The image is segemented in GrayMatter, WhiteMatter and CerebrospinalFluid.

.. code-block:: bash

    "MagneticResonanceImage":
        {
            "MaterialCoding": {
                "Unknown": 0,
                "GrayMatter": 1,
                "WhiteMatter": 2,
                "CerebrospinalFluid": 3
            },
            "Path": "./input_files/segmask.nii"
        },


For this study we choose a custom made monopolar electrode from Microprobes especially designed
for small rodents. 
Since the electrode is symetric the rotaion has no impact and is left to 0.
Further we want to implant the electrode crainial into the region of the subthalamic nucleus (STN)
relativ to the given MRI. This electrode has only one contact so we set it as active with a
voltage of 1 V.

.. code-block:: bash

    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "Rotation": 0,
            "Direction": [2.2314, 5.0661, 3.990],
            "Translation": [14.5884, -14.7434, -9.0634],
            "Contact_1": {"Active": true,
                          "Value": 1.0},
    ]


As a simplification we choose to set the skull of the animal as ground with a voltage of 0 V.
We can use this simplification since the distance from the electrode tip to the skull is relativly
large and the electric field is decreasing rapidly in the vincinity of the active contact.

.. code-block:: bash

    "BrainSurface":
        {
            "Active": true,
            "Value": 0
        },


As optional parameter we introduce a region of interest around our implantation centre. By setting the
size of the region we define an area in which the simulation will have a higher resolution.

.. code-block:: bash

    "RegionOfInterest":
        {
            "Active": true,
            "Center": [14.937, -13.613, -5.123],
            "Shape": [80, 80, 80]
        },

For the explanation of further parameters you can ferfer to the general introduction of the input files.

.. note::

    Format of the input needs to be updated!

An example of a full input directorie can be find below.

.. code-block:: bash

    {
    "BrainSurface":
        {
            "Active": false,
            "Value": 0
        },
    "DiffusionTensorImage":
        {
            "Path": ""
        },
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirected",
            "Rotation": 6.412,
            "Direction": [2.2314, 5.0661, 3.990],
            "Translation": [14.5884, -14.7434, -9.0634],
            "Contact_1": {"Active": true,
                          "Value": 1.0},
            "Contact_2": {"Active": false,
                          "Value": 0.0},
            "Contact_3": {"Active": false,
                          "Value": 0.0},
            "Contact_4": {"Active": false,
                          "Value": 0.0},
            "Contact_5": {"Active": false,
                          "Value": 0.0},
            "Contact_6": {"Active": false,
                          "Value": 0.0},
            "Contact_7": {"Active": false,
                          "Value": 0.0},
            "Contact_8": {"Active": true,
                          "Value": 0.0},

            "Contacts": {
                "Active": [true,
                           false,
                           false,
                           false,
                           false,
                           false,
                           false,
                           true],
                "Value": [1.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0]
            }
        }
    ],
    "FEMMode": "QS",
    "MagneticResonanceImage":
        {
            "MaterialCoding": {
                "Unknown": 0,
                "GrayMatter": 1,
                "WhiteMatter": 2,
                "CerebrospinalFluid": 3
            },
            "Path": "./input_files/segmask.nii"
        },
    "Mesh": {
        "LoadMesh": false,
        "LoadPath": "./input_files/mesh.vol",
        "MeshElementOrder": 2,
        "SavePath": ""
    },
    "MeshElementOrder": 2,
    "OutputPath": "test_result",
    "RegionOfInterest":
        {
            "Active": true,
            "Center": [14.937, -13.613, -5.123],
            "Shape": [80, 80, 80]
        },
    "SpectrumMode": "NoTruncation",
    "StimulationSignal":
        {
            "Type": "Rectangle",
            "Frequency": 130.0,
            "PulseWidthMicroSeconds": 60.0,
            "TopWidthMicroSeconds": 0.0
        }
