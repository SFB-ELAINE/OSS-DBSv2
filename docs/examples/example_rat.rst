Rodent Studies
==============

This example demonstrates how to perform a DBS simulation in a rodent brain using OSS-DBSv2.

MRI and DTI data are openly available, for example from the Johnson et al.\ rat brain atlas:

- **Publication**: Johnson, G. Allan, et al. *A multicontrast MR atlas of the Wistar rat brain.* NeuroImage, 242 (2021), 118470.
  DOI: https://doi.org/10.1016/j.neuroimage.2021.118470
- **Original Dataset**: https://civmvoxport.vm.duke.edu/voxbase/studyhome.php?studyid=754/
- **License**: https://creativecommons.org/licenses/by-nc-sa/3.0/
- **Author/Creator**: Johnson, G. Allan, et al.

After downloading the dataset, the segmented MRI (``segmask.nii.gz``) and the normalized DTI can be used in OSS-DBSv2 as follows:

.. code-block:: json

    "MaterialDistribution": {
      "MRIPath": "./input_files/segmask.nii.gz",
      "MRIMapping": {
        "Unknown": 0,
        "CSF": 3,
        "White matter": 2,
        "Gray matter": 1,
        "Blood": 4
      },
      "DiffusionTensorActive": false,
      "DTIPath": ""
    },

Rodent Electrodes
-----------------

In addition to clinical DBS electrodes, OSS-DBSv2 includes dedicated electrode models for small rodent studies.
A complete list is available in the :ref:`electrode documentation <electrodes>`.

In this example, a SNEX100 electrode is implanted cranially into the region of the subthalamic nucleus (STN).
The electrode has two contacts, so either bipolar or monopolar stimulation can be used.
Below a configuration for monopolar stimulation with 1 V at contact 1 is shown.
As simplification, the outer boundary of the brain region is treated as ground.
Since most rodent electrodes do not feature directional contacts, no orientation needs to be specified.

.. code-block:: json

    "Electrodes": [
      {
        "Name": "MicroProbesSNEX100",
        "Rotation[Degrees]": 0,
        "Direction": {"x[mm]": 2.23, "y[mm]": 5.07, "z[mm]": 3.99},
        "TipPosition": {"x[mm]": 14.59, "y[mm]": -14.74, "z[mm]": -9.06},
        "Contacts": [
            {
              "Contact_ID": 1,
              "Active": true,
              "Voltage[V]": 1.0,
              "Floating": false
            },
                            {
                    "Contact_ID": 2,
                    "Active": false,
                    "Voltage[V]": 0.0,
                    "Floating": true,
                }

        ]
      }
    ]
    "Surfaces": [
        {
            "Name": "BrainSurface",
            "Active": true,
            "Voltage[V]": 0.0,
            "Floating": false
        }
    ],

Estimate Stimulation Volume
---------------------------

To the stimulation volume a point model around the electrode tip is created.
Therefore, the ``Lattice`` option is activated and the location and number of points specified.

.. code-block:: json

  "PointModel": {
    "Pathway": {
      "Active": false,
      "FileName": ""
    },
    "Lattice": {
      "Active": true,
      "Center":  {"x[mm]": 14.59, "y[mm]": -14.74, "z[mm]": -9.06},
      "Shape": {"x": 20, "y": 20, "z": 20},
      "Direction": {"x[mm]": 2.23, "y[mm]": 5.07, "z[mm]": 3.99},
      "PointDistance[mm]": 0.5
    }
  },

The results are stored in the specified ``OutputPath`` directory after running the simulation.
