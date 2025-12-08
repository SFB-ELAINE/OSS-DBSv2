
Using OSS-DBSv2 with Lead-DBS
=============================

Overview
--------

Lead-DBS is the main entry point for most users of OSS-DBSv2.
It handles image preprocessing, atlas registration, and electrode reconstruction via its graphical user interface.
OSS-DBSv2 acts as a backend field solver and activation model that Lead-DBS calls in the background.

For an average Lead-DBS user, no separate OSS-DBSv2 setup is required.
The toolbox is installed automatically into a dedicated conda environment the first time you select OSS-DBSv2 as the electric field model in the Lead-DBS GUI.

Workflow via the Lead-DBS GUI
-----------------------------

A typical workflow for human DBS projects is:

1. Run the Lead-DBS pipeline for your patient:
   coregistration, normalization (optional), and electrode reconstruction.
2. Open the electrode scene and stimulation module in the Lead-DBS GUI.
3. Select **OSS-DBSv2** as the model for electric field computation.
4. Run the calculation for the stimulation volume or pathway activation.

Lead-DBS will then

- create the required input files in the patient’s ``stimulations`` folder,
- install or update the OSS-DBSv2 conda environment if needed, and
- call OSS-DBSv2 in the background to perform the simulation.

The resulting outputs are visualized directly in the electrode scene and stored in the patient’s folder.

Using OSS-DBSv2 directly
------------------------

For more advanced use cases (e.g., parameter sweeps, optimization, uncertainty quantification), you can use the Lead-DBS output as input to standalone OSS-DBSv2 runs.

Lead-DBS writes a stimulation-specific ``stimparameters.mat`` file into the patient’s ``stimulations`` subfolder.
This file contains the complete stimulation configuration in a format that can be converted into OSS-DBSv2 settings with:

.. code-block:: bash

   leaddbs2ossdbs path/to/stimparameters.mat

This command generates an OSS-DBSv2 configuration file in JSON format that can be used as input for simulations, as described in the tutorial.

Lead-DBS still provides all image and electrode information, while OSS-DBSv2 is controlled explicitly from the command line or from scripts for large-scale or customized simulations.

