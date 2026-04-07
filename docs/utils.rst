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
