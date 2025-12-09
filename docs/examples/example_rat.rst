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

In this example, the electrode is implanted cranially into the region of the subthalamic nucleus (STN).
This electrode has a single contact, which we activate at 1 V.
Since most rodent electrodes do not feature directional contacts, no orientation needs to be specified.

.. code-block:: json

    "Electrodes": [
      {
        "Name": "MicroProbesSNEX100",
        "Rotation": 0,
        "Direction": [2.2314, 5.0661, 3.990],
        "Translation": [14.5884, -14.7434, -9.0634],
        "Contact_1": {
          "Active": true,
          "Value": 1.0
        }
      }
    ]

Boundary Condition (Ground)
---------------------------

As a simplification, we treat the skull of the animal as ground (0 V).
This is acceptable because the distance between the electrode tip and skull is relatively large, and the electric field decays rapidly with distance.

.. code-block:: json

    "BrainSurface": {
      "Active": true,
      "Value": 0
    },
