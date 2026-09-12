Utilities
=========

This section groups smaller helper functionality used across the package. These
utilities are mostly relevant for developers, advanced scripting workflows, and
debugging.

Type checking
-------------

OSS-DBSv2 validates settings before running larger simulations. This reduces the
chance that expensive solves fail only after geometry construction or meshing.

Other utility modules support tasks such as:

- field post-processing
- material lookup and settings handling
- NIfTI image interaction
- VTK export

For new users, these helpers usually stay in the background. They become more
important when extending the package or building custom automation around it.

Related pages
-------------

- :doc:`python_api` — scripting API for driving OSS-DBSv2 from Python
- :doc:`input_settings` — JSON settings reference (uses ``Settings`` defaults)
- :doc:`volume_conductor_model` — the FEM stack that uses these utilities

API reference
-------------

Settings and defaults
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.utils.settings
    :members:
    :undoc-members:
    :show-inheritance:

Type checking
^^^^^^^^^^^^^

.. automodule:: ossdbs.utils.type_check
    :members:
    :undoc-members:
    :show-inheritance:

NIfTI image handling
^^^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.utils.nifti1image
    :members:
    :undoc-members:
    :show-inheritance:

VTK export
^^^^^^^^^^

.. automodule:: ossdbs.utils.vtk_export
    :members:
    :undoc-members:
    :show-inheritance:

Field computation
^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.utils.field_computation
    :members:
    :undoc-members:
    :show-inheritance:

VTA collapse
^^^^^^^^^^^^

.. automodule:: ossdbs.utils.collapse_vta
    :members:
    :undoc-members:
    :show-inheritance:

Materials
^^^^^^^^^

.. automodule:: ossdbs.utils.materials
    :members:
    :undoc-members:
    :show-inheritance:
