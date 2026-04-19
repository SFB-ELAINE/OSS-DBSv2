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
- A multi-contact example with ``C1: -3 mA, C2: 1 mA, C3: 1 mA, C4: 1 mA``.

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

This case demonstrates the **StimSets** batch workflow, which computes
unit solutions for each electrode contact and then superimposes them for
an arbitrary number of stimulation protocols.

**Setup:**

- A Medtronic 3387 electrode (4 contacts) is placed in an ellipsoidal
  brain region.
- The ``Surfaces`` block sets ``BrainSurface`` as a current-controlled
  ground (``Current[A]: -1``).
- Stimulation is current-controlled with a 130 Hz rectangular pulse.
- ``StimSets`` is enabled and points to ``Current_protocol.csv``.

**How it works:**

1. The mesh is generated once and saved (with material bisection applied).
2. For each non-ground contact, the code creates a unit-current solve:
   only that contact is active (1 A), all others are floating, and the
   ground carries −1 A.
3. HP refinement is re-applied on the loaded mesh for each unit solve.
4. The unit solutions are stored in per-contact output directories
   (``E1C1``, ``E1C2``, …).

**Pathway activation (PAM):**

After the unit solutions are computed, pathway activation is run
separately via ``run_pathway_activation``:

.. code-block:: bash

   run_pathway_activation input_pathway.json

This loads all unit solutions, reads each row of ``Current_protocol.csv``
as a scaling vector, superimposes the unit solutions accordingly, and
runs the NEURON model for each protocol. See :ref:`stimsets-workflow` in
the settings reference for full details.

**CSV format:**

``Current_protocol.csv`` is a CSV file with a header row naming the
contacts and one row per stimulation protocol. Values are in
**milliamps** (converted to Amperes internally). NaN entries are treated
as zero current.

.. code-block:: text

   Contact0,Contact1,Contact2,Contact3
   -1.954,0.0,-0.983,0.0
   -3.811,-0.933,-2.437,0.784
   1.113,-1.445,3.622,-3.436


Case 10: Floating with Surface Impedance
-----------------------------------------

This case demonstrates the ``VolumeConductorFloatingImpedance`` formulation
using a PINS Medical L303 electrode with three contacts.

Three input files are provided:

- ``input_floating.json``: plain floating contacts (no surface impedance).
  Contact 1 at 1 V, contact 3 at ground, contact 2 floating.
- ``input_floating_impedance.json``: same contact activation but all floating
  contacts carry a resistive surface impedance (``R = 1 kOhm``).
- ``input_floating_impedance_noground.json``: all three contacts are floating
  with surface impedance and no active ground. Currents are prescribed
  (+1 mA / 0 / -1 mA) in current-controlled mode. A Lagrange multiplier
  constrains the sum of floating potentials to zero.

This case is useful for verifying that:

- floating contacts with surface impedance produce physically plausible results
- the no-ground current-controlled configuration is well-posed
- impedance estimates include the surface-impedance dissipation

See :ref:`surface-impedance` and :ref:`stimulation-modes` for background.


Case 11: Impedance Validation (Homogeneous)
--------------------------------------------

This case validates impedance computation against analytical or reference
values in a homogeneous tissue medium using a Boston Scientific Vercise
electrode.

Six input files explore different combinations:

- ``input_homogeneous.json``: baseline without interface impedance
- ``input_no_interface.json``: heterogeneous MRI-based tissue, no interface
- ``input_homogeneous_lempka2009.json``: homogeneous tissue with Lempka (2009)
  CPE interface parameters (``CPE_dl``, ``dl_k = 1.5e6``, ``dl_alpha = 0.8``)
- ``input_interface_lempka2009.json``: MRI-based tissue with Lempka interface
- ``input_surface_impedance_lempka2009.json``: surface impedance (Robin BC)
  formulation with Lempka parameters on active contacts
- ``input_homogeneous_surface_impedance_lempka2009.json``: homogeneous tissue
  with surface impedance

These configurations allow comparing the lumped-element interface model with
the surface impedance Robin BC approach, and verifying that both give
consistent impedance spectra.


