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
- ``Shape``: for example ``Sphere`` or another supported geometry type

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

Important solver settings include:

- ``Type`` such as ``CG``, ``GMRES``, or ``Direct``
- ``Preconditioner`` such as ``bddc`` or ``local``
- ``MaximumSteps`` and ``Precision``

For most new users, the example defaults are a good starting point. Solver and
mesh tuning usually becomes important only for large studies or demanding
anisotropic models.

Stimulation signal
------------------

``StimulationSignal`` defines how the excitation is represented. Common fields
include:

- ``CurrentControlled`` to switch between current- and voltage-controlled setups
- ``Type`` for the signal class
- frequency or pulse-shape parameters depending on the chosen model

Different workflows use different subsets of these settings. When importing
from Lead-DBS, several values are generated automatically.

Point models and post-processing
--------------------------------

``PointModel`` controls field evaluation away from the FEM mesh:

- ``Lattice`` evaluates the field on a regular grid and can be used for VTA-like
  analyses
- ``Pathway`` activates pathway-based analysis when tract data are available

These options are often disabled for a first validation run and enabled later
once the core volume conductor setup is working as expected.

Output files
------------

``OutputPath`` defines where OSS-DBSv2 writes results. In addition to numerical
outputs, the software writes log files and status flags that are especially
useful in automated or Lead-DBS-driven workflows.

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
