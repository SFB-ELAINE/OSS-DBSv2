Volume conductor model
======================

The volume conductor model is the numerical core of OSS-DBSv2. It combines the
geometry, tissue properties, stimulation signal, and solver configuration to
compute electric potentials and fields in the simulated domain.

Role in the workflow
--------------------

Most standalone and Lead-DBS-driven runs eventually reach the same core steps:

1. build the geometry and mesh
2. assign conductivity information from MRI and optional DTI data
3. define boundary conditions at electrode contacts and outer surfaces
4. solve the resulting finite-element system
5. export fields, impedance estimates, and pointwise analyses

Floating and non-floating contacts
----------------------------------

OSS-DBSv2 distinguishes between several conductor variants:

- ``VolumeConductorNonFloating`` for standard active and passive contact setups
- ``VolumeConductorFloating`` when floating contacts are present
- ``VolumeConductorFloatingImpedance`` when floating contacts also carry surface
  impedance information

These variants share the same general workflow but differ in how contact
boundary conditions are enforced.

Mesh and linear algebra
-----------------------

The volume conductor stack also includes:

- mesh generation and mesh reuse utilities
- conductivity coefficient functions derived from the material model
- iterative and direct solvers
- several preconditioner strategies for larger problems

For most users, the default settings are sufficient to begin with. Solver and
preconditioner tuning becomes more relevant for large anisotropic models, high
resolution meshes, or systematic parameter studies.

API reference
-------------

Volume conductor model
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.volume_conductor.volume_conductor_model
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.nonfloating
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.floating
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.fem.volume_conductor.floating_impedance
    :members:
    :undoc-members:
    :show-inheritance:

Material parameters
^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.volume_conductor.conductivity
    :members:
    :undoc-members:
    :show-inheritance:


Mesh
^^^^

.. automodule:: ossdbs.fem.mesh
    :members:
    :undoc-members:
    :show-inheritance:

Solver
^^^^^^

.. automodule:: ossdbs.fem.solver
    :members:
    :undoc-members:
    :show-inheritance:


Preconditioner
^^^^^^^^^^^^^^

.. automodule:: ossdbs.fem.preconditioner
    :members:
    :undoc-members:
    :show-inheritance:
