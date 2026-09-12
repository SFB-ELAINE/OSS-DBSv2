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
assignment (see :ref:`dti-anisotropy` below).

.. _dti-anisotropy:

Anisotropic conductivity from DTI
---------------------------------

Diffusion tensor imaging (DTI) data can be used to model the directional
dependence of tissue conductivity, which is particularly important in white
matter where myelinated fibre bundles cause strongly anisotropic conduction.

File format
^^^^^^^^^^^

The DTI input must be a 4-D NIfTI file (`.nii` or `.nii.gz`) with shape
`(x, y, z, 6)`. The last dimension stores the six unique components of the
symmetric 3×3 diffusion tensor in the order:

.. code-block:: text

   index 0: xx
   index 1: xy  (= yx)
   index 2: xz  (= zx)
   index 3: yy
   index 4: yz  (= zy)
   index 5: zz

The file must use the same spatial reference frame as the segmented MRI.
Spatial units in the NIfTI header are respected (meter, mm, or micron are
converted to mm automatically).

The DTI data needs to be preprocessed before importing.
Currently, this is supported in Lead-DBS where common preprocessing steps are implemented.
The Lead-DBS template data provides an example
(``IITmean_tensor_Norm_mapping.nii.gz``).

JSON configuration
^^^^^^^^^^^^^^^^^^

DTI is enabled in the ``MaterialDistribution`` block:

.. code-block:: json

   "MaterialDistribution": {
     "MRIPath": "segmask.nii.gz",
     "MRIMapping": {
       "Unknown": 0,
       "CSF": 1,
       "White matter": 2,
       "Gray matter": 3,
       "Blood": 4
     },
     "DiffusionTensorActive": true,
     "DTIPath": "dti_norm_mapping.nii.gz",
     "WMMasking": true
   }

- ``DiffusionTensorActive`` (default: ``false``) — set to ``true`` to load
  and use DTI data.
- ``DTIPath`` — path to the DTI NIfTI file.
- ``WMMasking`` (default: ``true``) — when ``true``, the DTI tensor is
  applied only in white-matter voxels; all other tissues receive an
  isotropic identity tensor. When ``false``, the DTI tensor is applied
  everywhere.

How conductivity is constructed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At each voxel and each frequency, the final conductivity is computed as:

.. math::

   \boldsymbol{\sigma}(\mathbf{x}, \omega) =
   \begin{cases}
     \mathbf{D}(\mathbf{x}) \cdot \sigma_{\text{scalar}}(\mathbf{x}, \omega)
       & \text{if white matter (or WMMasking is off)} \\
     \mathbf{I} \cdot \sigma_{\text{scalar}}(\mathbf{x}, \omega)
       & \text{otherwise}
   \end{cases}

where :math:`\mathbf{D}(\mathbf{x})` is the normalised DTI tensor,
:math:`\sigma_{\text{scalar}}` is the scalar conductivity from the
dielectric model (already converted to S/mm), and :math:`\mathbf{I}` is
the 3×3 identity tensor.

This means the DTI tensor acts as a directional scaling of the scalar
tissue conductivity. Encapsulation layers are always treated as isotropic,
even when DTI is active.

When DTI is enabled, the conductivity is returned as a tensor-valued
(3×3 matrix) coefficient function. All downstream FEM operations
(assembly, export) handle this automatically.

Practical guidance
^^^^^^^^^^^^^^^^^^

- White-matter masking (``WMMasking: true``) is recommended. Without it,
  the DTI tensor is applied in all tissues, which can produce unphysical
  anisotropy in CSF or gray matter.
- Verify the DTI and MRI images share the same spatial reference. Mis-
  aligned images will map the wrong tensor to each voxel.
- The ``examples/ImageAPI/example_dti_mask.py`` script demonstrates how
  to inspect masked vs. unmasked conductivity fields and export them for
  visualisation in ParaView.

Practical guidance
------------------

- Keep the MRI label mapping consistent with your segmentation pipeline.
- Use built-in defaults first before introducing custom parameter sets.
- When comparing studies, document both the dielectric model family and any
  custom parameters, since they can materially influence the results.

Related pages
-------------

- :doc:`input_settings` — ``DielectricModel`` and ``MaterialDistribution`` settings
- :doc:`brain_geometry` — how tissue regions are defined in the geometry
- :doc:`volume_conductor_model` — how conductivity enters the FEM formulation

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
