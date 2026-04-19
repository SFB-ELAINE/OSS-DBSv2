Axon models
===========

OSS-DBSv2 supports pathway-based analyses in which extracellular fields are
sampled along trajectories and then coupled to axon models. This extends the
workflow beyond field simulation toward estimates of neural pathway activation.

Workflow overview
-----------------

A typical pathway-based workflow consists of:

1. generating or importing pathway trajectories
2. evaluating the extracellular solution along those pathways
3. mapping the sampled signal to a neuron model
4. computing activation-related outputs

This part of the software is more specialized than the basic volume conductor
workflow, but it becomes essential when the scientific question focuses on
which pathways are activated rather than only how far the field spreads.

Implemented model families
--------------------------

The codebase includes morphology handling for models such as MRG2002 and
McNeal1976, together with a neuron-simulation abstraction for running pathway
activation analyses.

When to use this part of the toolbox
------------------------------------

- use it when tract-based activation is a central study outcome
- skip it for a first validation run focused only on electric fields or basic
  VTA-like analyses
- treat it as a separate workflow layer that builds on a correct and stable
  volume conductor setup

Related pages
-------------

- :doc:`point_analysis` — lattice and pathway evaluation that feeds into
  axon models
- :doc:`volume_conductor_model` — the upstream FEM solve
- :doc:`input_settings` — ``PointModel.Pathway`` configuration

API reference
-------------


Axon models
^^^^^^^^^^^

.. automodule:: ossdbs.axon_processing.axon_models
    :members:
    :undoc-members:
    :show-inheritance:

Neuron simulator
^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.axon_processing.neuron_model
    :members:
    :undoc-members:
    :show-inheritance:
