Brain geometry
==============

OSS-DBSv2 distinguishes between a brain-only geometry and the final simulation
geometry. The brain geometry defines the simulation domain, while the full model
geometry combines that domain with implanted electrodes, contact surfaces, and
optional encapsulation layers.

Two geometry levels
-------------------

- ``BrainGeometry`` describes the outer brain domain used for meshing and field
  computation.
- ``ModelGeometry`` combines the brain domain with one or more electrode models
  and assigns the contact and material regions needed for the FEM solve.

In practical terms, most users configure geometry indirectly through the input
JSON. The code then constructs the necessary CAD objects automatically.

Brain geometry definition
-------------------------

The brain domain can be defined in two main ways:

- as a simple analytic shape such as a sphere, ellipsoid, or box
- from imaging-derived information, for example by using the MRI bounding box
  and affine transformation

This level is mainly responsible for setting the spatial extent of the
simulation and naming the outer surface as ``BrainSurface``.

Electrode integration
---------------------

Electrodes are created separately and then inserted into the brain domain. When
the final ``ModelGeometry`` is assembled, OSS-DBSv2:

- places the electrode CAD models according to position, direction, and rotation
- creates contact surfaces with stable names such as ``E1C1``
- adds encapsulation layers when requested
- prepares surface and volume names used later for materials, boundary
  conditions, and mesh refinement

This step is central because it connects the anatomical model to the numerical
problem that is solved by the volume conductor model.

Practical considerations
------------------------

- Choose a brain region that is large enough to avoid boundary effects but not
  so large that meshing becomes unnecessarily expensive.
- Start with standard electrode models whenever possible before introducing
  custom geometries.
- If geometry construction fails, reduce complexity first: smaller domain,
  simpler shape, or no encapsulation layer.
- Use the example cases to understand how contact numbering maps to the final
  geometry.

API reference
-------------

Brain geometry
^^^^^^^^^^^^^^

.. automodule:: ossdbs.model_geometry.brain_geometry
    :members:
    :undoc-members:
    :show-inheritance:

Model geometry
^^^^^^^^^^^^^^

.. automodule:: ossdbs.model_geometry.model_geometry
    :members:
    :undoc-members:
    :show-inheritance:


Helper classes
^^^^^^^^^^^^^^

.. automodule:: ossdbs.model_geometry.bounding_box
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.model_geometry.contacts
    :members:
    :undoc-members:
    :show-inheritance:
