OSS-DBSv2 Input Test Cases
==========================

This section documents the example input test cases provided with OSS-DBSv2.
They demonstrate different modeling features, parameter settings, electrode configurations, and output options.

Running a Test Case
-------------------

Each test case consists of an input file in JSON format stored in the ``input_test_cases`` directory.
To run a test case:

.. code-block:: bash

    ossdbs path/to/input_file.json

During execution, diagnostic output is printed to the console.
Upon completion, results are written to an automatically created output folder located (by default) next to the JSON file.

Depending on the configuration, the results folder may contain:

- Material and conductivity distribution (``.vtk``)
- Electric potential and electric field distribution (``.vtk``)
- Probed electric potential and field (``.h5``)
- Volume of tissue activated (VTA) (``.nii``)
- Electrode geometry (``.vtk``)
- Generated mesh (``.vol.gz``)
- Estimated impedance spectrum (``.csv``)


Case 1: Brain Tissue Properties
-------------------------------

This case demonstrates the use of homogeneous vs.\ inhomogeneous and anisotropic tissues.

Two input files are provided:

- First setup: uses a homogeneous MRI segmentation, no diffusion tensors.
- Second setup: uses a real MRI scan of a human brain.

A Boston Scientific Vercise electrode is placed in the STN.  
The lowest contact is driven with 1 V, the second remains at ground potential.

The dielectric properties follow the ColeCole4 model at 10 kHz, using the quasi-static (QS) approximation.  
All results are written into the corresponding ``results`` folder.


Case 2: Custom Parameters
-------------------------

This case demonstrates how to modify electrode geometry and dielectric properties.

- The directed Boston Scientific Vercise electrode is placed in the STN.
- A 0.2 mm encapsulation layer is added around the shaft (default: gray matter).
- In ``input_custom_electrode.json``, the electrode is altered by adjusting contact lengths.
- A second input file shows how to override dielectric material properties with constants.

This example highlights how flexible user-defined parameterization can be.


Case 3: Case Grounding
----------------------

This case uses a custom rodent monopolar electrode (MicroProbes) implanted in the STN, using the same MRI as in Case 1.

- The electrode drives 1 V stimulation.
- The outer boundary of the brain region is treated as ground.

For custom geometries, any part of the boundary or specified regions can be treated as ground analogously.


Case 4: Current-Controlled Stimulation
--------------------------------------

This case demonstrates stimulation in **current-controlled** mode.

- An Abbott St. Jude directed electrode is used.
- Two example configurations are provided:
  - A bipolar configuration with ±0.1 mA.
  - A multi-contact example:
    ``C1: –3 mA, C2: 1 mA, C3: 1 mA, C4: 1 mA``

Multiple current-driven contacts can be stimulated simultaneously.


Case 5: Stimulation Signals
---------------------------

This case exemplifies defining and solving **time-domain stimulation signals**.

- A Medtronic SenSight electrode is used.
- The rectangular signal has a base frequency of 130 Hz and a pulse width of 60 µs.
- The signal is transformed into the frequency domain, solved, and converted back to time-domain response.
- The octave-band method is used to reduce the number of sampled frequencies.


Case 6: Floating Contacts
-------------------------

Floating conductors are demonstrated using the PINS Medical L303 electrode.

- Contact 1: 1 V stimulation.
- Contact 3: Ground.
- Contact 2: Floating (electrically passive).

Floating contacts influence the field solution without being explicitly driven.


Case 7: Volume of Tissue Activated (VTA)
----------------------------------------

This test case estimates VTA using point evaluation around the active contact of a DIXI SEEG electrode.

- A uniform grid of evaluation points is defined.
- The electric field is probed and thresholded.
- The resulting VTA is stored in NIfTI format.

Optionally, the electrode volume can be removed from the VTA (collapsing it inside the trajectory).


Case 8: Pathway Activation Modelling (PAM)
------------------------------------------

This example shows pathway-based evaluation of electric fields.

- A Medtronic 3387 electrode is used in the STN.
- Field values are sampled along axon trajectories defined in an ``.h5`` file.
- The ``.h5`` file must contain groups (pathways) containing datasets (axons):
  ``axon0``, ``axon1``, ``axon2``, etc.
- Each dataset contains a 3×N array of spatial coordinates.
- Results are stored in ``.h5`` format using the same structure.


Case 9: StimSets
----------------

WIP – will demonstrate batch simulations of the same electrode configuration across varying stimulation settings.

