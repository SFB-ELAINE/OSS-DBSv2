Point analysis
==============

Point analysis evaluates the FEM solution at selected locations after the volume
conductor problem has been solved. This is the bridge between the raw field
solution and downstream analyses such as VTA-like thresholding or pathway-based
activation studies.

Available point models
----------------------

OSS-DBSv2 includes multiple ways to sample the field:

- ``Lattice`` for regular point grids around the electrode
- ``VoxelLattice`` for image-like outputs aligned with voxel space
- ``Pathway`` for points that lie on fiber trajectories or axon pathways

These models are especially important when the primary output is not just the
potential field itself, but a derived spatial analysis.

Typical use cases
-----------------

- estimate VTA-like regions from electric field thresholds
- export pointwise potentials or field magnitudes for post-processing
- evaluate pathway trajectories for later axon activation modeling

Practical guidance
------------------

- Disable point models during the very first setup run if you want to validate
  the core FEM workflow first.
- Enable lattices when exploring field spread around the electrode.
- Enable pathway models when tract data and axon activation workflows are part
  of the study design.


API reference
-------------

.. automodule:: ossdbs.point_analysis.point_model
    :members:
    :undoc-members:
    :show-inheritance:
    :noindex:

.. automodule:: ossdbs.point_analysis.lattice
    :members:
    :undoc-members:
    :show-inheritance:
    :noindex:

.. automodule:: ossdbs.point_analysis.pathway
    :members:
    :undoc-members:
    :show-inheritance:
    :noindex:

.. automodule:: ossdbs.point_analysis.voxel_lattice
    :members:
    :undoc-members:
    :show-inheritance:
    :noindex:

.. automodule:: ossdbs.point_analysis.time_results
    :members:
    :undoc-members:
    :noindex:
