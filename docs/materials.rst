Dielectric properties
=====================

The dielectric model determines how tissue conductivity and permittivity are
assigned during the simulation. In OSS-DBSv2, tissue labels from the MRI are
mapped to material names such as gray matter, white matter, CSF, or blood, and
the selected dielectric model converts those materials into frequency-dependent
electrical properties.

Why this matters
----------------

The electric field depends strongly on tissue properties. This page covers the
models used to translate anatomy into conductivity values for the volume
conductor solve.

Available model families
------------------------

OSS-DBSv2 currently provides:

- ``ColeCole4`` for frequency-dependent tissue behavior using a multi-pole
  Cole-Cole formulation
- ``ColeCole3`` as a related reduced model
- ``Constant`` for fixed conductivity and permittivity values

For many subject-specific workflows, a Cole-Cole model is the natural starting
point. The constant model is often useful for controlled comparisons, simplified
test cases, or custom homogeneous setups.

How material assignment works
-----------------------------

The workflow is typically:

1. provide a segmented MRI volume
2. map voxel labels to tissue names in ``MRIMapping``
3. select a dielectric model in ``DielectricModel``
4. optionally provide custom parameter overrides

If DTI is active, anisotropy can be incorporated on top of the scalar material
assignment.

Practical guidance
------------------

- Keep the MRI label mapping consistent with your segmentation pipeline.
- Use built-in defaults first before introducing custom parameter sets.
- When comparing studies, document both the dielectric model family and any
  custom parameters, since they can materially influence the results.

API reference
-------------

.. automodule:: ossdbs.dielectric_model.dielectric_model
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.dielectric_model.colecole4
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.dielectric_model.colecole3
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.dielectric_model.constant
    :members:
    :undoc-members:
    :show-inheritance:
