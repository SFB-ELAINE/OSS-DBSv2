Input and Settings Reference
============================

OSS-DBSv2 simulations are configured through a JSON input file. This page gives
an overview of the most important top-level sections so external users can map
their modeling question to the required settings.

For a runnable minimal example, see
`input_files/inputTest.json <https://github.com/SFB-ELAINE/OSS-DBSv2/blob/main/input_files/inputTest.json>`_.
Use that file together with :doc:`tutorial` if you want to compare the
reference explanations here with a working end-to-end example.

Top-level structure
-------------------

The most common sections are:

- ``BrainRegion``: size, shape, and position of the simulated brain domain
- ``Electrodes``: electrode models, orientation, contact activation, and
  encapsulation layer settings
- ``MaterialDistribution``: MRI and optional DTI inputs plus tissue label mapping
- ``DielectricModel``: electrical material model used for each tissue class
- ``Mesh``: mesh generation and mesh-size controls
- ``StimulationSignal``: current- or voltage-controlled excitation settings
- ``Solver``: linear solver and preconditioner configuration
- ``PointModel``: lattice- or pathway-based post-processing
- ``OutputPath``: destination for logs and generated result files

Not every workflow needs every section. A first standalone run usually depends
most on ``BrainRegion``, ``Electrodes``, ``MaterialDistribution``,
``StimulationSignal``, and ``OutputPath``.

Brain region
------------

``BrainRegion`` defines the simulation domain around the electrode. A typical
section contains:

- ``Center``: the center of the domain in millimeters
- ``Dimension``: the extent in x, y, and z direction
- ``Shape``: for example ``Sphere`` or another supported geometry type like ``Box`` or ``Ellipsoid`` 

For standalone usage, this is often the first section to adapt when moving from
an example dataset to a new subject or experiment.

Electrodes
----------

``Electrodes`` is a list. Each entry describes one implanted electrode with:

- ``Name``: the electrode model to instantiate
- ``TipPosition``: electrode tip location in millimeters
- ``Direction``: implantation axis
- ``Rotation[Degrees]``: rotation around the electrode axis
- ``Contacts``: activity and boundary conditions for each contact
- ``EncapsulationLayer``: optional peri-electrode layer properties

Each contact can specify whether it is active, floating, current-controlled, or
voltage-controlled. Contact-level mesh settings can also be added when local
refinement is needed.

Material distribution
---------------------

``MaterialDistribution`` connects imaging data to tissue classes:

- ``MRIPath`` points to a segmented MRI volume
- ``MRIMapping`` translates integer labels in the MRI to OSS-DBSv2 tissue names
- ``DiffusionTensorActive`` enables anisotropy from DTI data
- ``DTIPath`` points to the diffusion tensor image when used

This section is central for subject-specific modeling. If ``DiffusionTensorActive``
is ``false``, the simulation uses the scalar tissue properties only.

Dielectric model
----------------

``DielectricModel`` selects how tissue conductivity and permittivity are modeled.
The code currently exposes models such as ``ColeCole4``, ``ColeCole3``, and
``Constant``.

Typical usage:

- use a Cole-Cole model for frequency-dependent tissue properties
- use the constant model when a homogeneous or fixed-conductivity setup is needed
- provide ``CustomParameters`` when deviating from the built-in defaults

Mesh and solver
---------------

``Mesh`` and ``Solver`` control numerical performance and accuracy.

Important mesh settings include:

- ``LoadMesh`` / ``SaveMesh`` for mesh reuse
- ``MeshElementOrder`` for the geometric discretization
- ``MeshingHypothesis`` for coarse-to-fine global mesh sizing
- ``MeshSize`` for targeted local refinement
- ``HPRefinement`` for combined geometric and polynomial refinement near
  boundaries

See :ref:`mesh-refinement` for a full description of all mesh control options.

Important solver settings include:

- ``Type`` such as ``CG``, ``GMRES``, or ``Direct``
- ``Preconditioner`` such as ``bddc`` or ``local``
- ``MaximumSteps`` and ``Precision``

For most new users, the example defaults are a good starting point. Solver and
mesh tuning usually becomes important only for large studies or demanding
anisotropic models. See :ref:`solver-guidance` for detailed recommendations on
solver and preconditioner selection.

Stimulation signal
------------------

``StimulationSignal`` defines how the excitation is represented. Common fields
include:

- ``CurrentControlled`` to switch between current- and voltage-controlled setups
- ``Type`` for the signal class
- frequency or pulse-shape parameters depending on the chosen model

Different workflows use different subsets of these settings. When importing
from Lead-DBS, several values are generated automatically. See
:ref:`stimulation-modes` for a detailed description of how current- and
voltage-controlled stimulation interact with contact configurations.

Point models and post-processing
--------------------------------

``PointModel`` controls field evaluation away from the FEM mesh:

- ``Lattice`` evaluates the field on a regular grid and can be used for VTA-like
  analyses
- ``Pathway`` activates pathway-based analysis when tract data are available

These options are often disabled for a first validation run and enabled later
once the core volume conductor setup is working as expected.

Surfaces
--------

The optional ``Surfaces`` block sets boundary conditions on the outer brain
surface. This is primarily used to treat the brain boundary as a grounded
electrode in monopolar stimulation setups:

.. code-block:: json

   "Surfaces": [
     {
       "Name": "BrainSurface",
       "Active": true,
       "Voltage[V]": 0.0,
       "Floating": false
     }
   ]

Each entry uses the same fields as a contact definition (``Active``,
``Voltage[V]``, ``Floating``, ``Current[A]``). The ``Name`` must match a
boundary name in the geometry — the default outer boundary is called
``BrainSurface``. If ``Surfaces`` is omitted or empty, the outer boundary
is left as a natural (zero-flux) boundary condition.

.. _stimsets-workflow:

StimSets (batch stimulation)
----------------------------

The ``StimSets`` workflow enables efficient batch simulation of many
stimulation protocols on the same electrode geometry. Instead of running a
full mesh-generation and FEM solve for every protocol, StimSets exploits
linearity: it computes one **unit solution** per contact and then
superimposes them with protocol-specific scaling factors.

JSON configuration
^^^^^^^^^^^^^^^^^^

.. code-block:: json

   "StimSets": {
     "Active": true,
     "StimSetsFile": "Current_protocol.csv"
   }

- ``Active`` (default: ``false``) — enable the StimSets workflow.
- ``StimSetsFile`` — path to a CSV file containing stimulation protocols
  (see format below). Can also be ``null`` when using the Python API with
  a ``CurrentVector`` instead.

Requirements
^^^^^^^^^^^^

- Stimulation must be **current-controlled** (``CurrentControlled: true``).
- Exactly one contact or surface must serve as **ground** with
  ``Current[A]: -1``. This is typically the brain surface
  (``Surfaces`` block) or one active electrode contact.
- The preconditioner is automatically set to ``local`` (floating contacts
  require it).

CSV format
^^^^^^^^^^

The CSV file has a header row with contact names and one row per
stimulation protocol. Values are in **milliamps** (the code converts to
Amperes by multiplying by 1e-3). NaN values are treated as zero current.

.. code-block:: text

   Contact0,Contact1,Contact2,Contact3
   -1.954,0.0,-0.983,0.0
   -3.811,-0.933,-2.437,0.784

The number of columns must match the number of non-ground contacts.

Execution workflow
^^^^^^^^^^^^^^^^^^

The StimSets pipeline proceeds in two stages:

**Stage 1 — Unit solutions** (``ossdbs``):

1. A single mesh is generated and h-refined (material bisection), then
   saved to disk.
2. For each non-ground contact, the mesh is reloaded, HP refinement is
   re-applied, and a unit-current FEM solve is performed (1 A on that
   contact, all others floating, ground at −1 A).
3. Unit solutions are stored in per-contact directories named
   ``<OutputPath>E1C1``, ``<OutputPath>E1C2``, etc.

**Stage 2 — Pathway activation** (``run_pathway_activation``):

1. All unit solutions are loaded from the per-contact directories.
2. For each row in the CSV file, the unit solutions are scaled by the
   protocol currents and superimposed.
3. The resulting time-domain field is passed to the NEURON model for
   axon activation analysis.

Stage 2 is run separately after stage 1:

.. code-block:: bash

   run_pathway_activation input.json

This separation allows re-running the PAM stage with different protocol
files or scaling factors without repeating the expensive FEM solves.

Additional top-level settings
-----------------------------

Several top-level keys control the FEM formulation and output behaviour:

- ``EQSMode`` (default: ``false``) — enable electro-quasi-static mode. When
  ``true``, the FEM solve uses complex-valued spaces to account for
  capacitive tissue effects at non-zero frequencies. When ``false``, only real
  conductivity is used (quasi-static approximation).
- ``FEMOrder`` (default: 2) — polynomial order of the finite-element space.
  Higher orders improve accuracy at the cost of more degrees of freedom.
  This is independent of ``MeshElementOrder``, which controls the geometric
  mesh curvature.
- ``ComputeImpedance`` (default: ``false``) — compute the complex impedance
  at each frequency and write ``impedance.csv``.
- ``ComputeCurrents`` (default: ``false``) — estimate contact currents at
  each frequency by integrating the normal current density over contact
  surfaces.
- ``ExportVTK`` (default: ``false``) — export potential, E-field,
  conductivity, and material distributions as VTU files for ParaView.
- ``ExportElectrode`` (default: ``false``) — export electrode geometry as VTU
  and Netgen mesh files.
- ``ExportFrequency`` (default: ``null``) — frequency at which VTK export is
  performed. If ``null``, the median frequency is used.

Output files
------------

``OutputPath`` defines where OSS-DBSv2 writes results. The following files may
be produced depending on the configuration:

**Always written:**

- ``ossdbs.log`` — full log of the simulation run
- ``VCM_report.json`` — degrees of freedom, element count, and solver timings

**Controlled by JSON flags:**

- ``ComputeImpedance: true`` produces ``impedance.csv`` with columns
  ``freq``, ``real``, ``imag``. For current-controlled multicontact setups,
  ``admittance_matrix.csv`` and ``impedance_matrix.csv`` are also written
  (flat format: ``freq``, ``row``, ``col``, ``real``, ``imag``).
- ``ExportVTK: true`` produces VTU files for ParaView visualisation:
  ``potential.vtu``, ``E-field.vtu``, ``conductivity.vtu``, and
  ``material.vtu``.
- ``ExportElectrode: true`` produces ``electrode_N.vtu`` and
  ``electrode_N.vol.gz`` for each electrode.

**Produced by specific workflows:**

- **Floating contacts**: ``floating_potentials.csv`` (frequency-domain floating
  voltages) and ``floating_in_time.csv`` (reconstructed time-domain values).
- **Time-domain signal**: ``stimulation_in_time.csv`` (reconstructed
  stimulation waveform at contact surfaces) and ``time_domain_signal.pdf``
  (plot of the reconstructed signal).
- **Point model (Lattice)**: ``oss_lattice_*.h5`` with potential and field
  values on the evaluation grid.
- **Point model (Pathway)**: ``oss_pathway_*.h5`` with field values sampled
  along axon trajectories, plus activation results when NEURON is available.
- **VTA**: ``oss_vta.nii`` (NIfTI volume of tissue activated).
- ``SaveMesh: true``: ``mesh.vol.gz`` for later reuse.

When adapting examples, it is often helpful to choose an explicit output folder
per subject or experiment so results from different runs do not get mixed.

Practical advice
----------------

- Start from a working example instead of writing a JSON file from scratch.
- Change one section at a time when adapting to a new dataset.
- Keep paths explicit when running batch jobs or working outside the example
  directories.
- Validate the electrode definition and MRI label mapping early, since these
  settings strongly influence downstream results.

Related pages
-------------

- :doc:`tutorial`
- :doc:`lead_dbs`
- :doc:`examples`
- :doc:`brain_geometry`
- :doc:`materials`
- :doc:`stimulation_signals`
- :doc:`volume_conductor_model`
