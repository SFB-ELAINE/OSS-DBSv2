.. OSS-DBS documentation master file, created by sphinx-quickstart on Fri Oct 21 11:27:48 2022.

Welcome to OSS-DBSv2's documentation!
=====================================

.. note::
   This page is still under construction.

Overview
--------

Deep brain stimulation (DBS) is a widely used treatment for several motor and non-motor disorders.
Because the exact mechanisms of action are still not fully understood, computational models are used to predict and optimize treatment outcomes.

OSS-DBSv2 is an open-source toolbox for highly automated DBS modeling in humans and animal models.
Users can provide their own magnetic resonance imaging (MRI) and diffusion-tensor imaging (DTI) data or use the included example datasets.
A set of predefined DBS electrodes covers the most common commercial designs, and additional custom electrodes can be added easily.
Stimulation parameters such as active contacts, amplitudes, and pulse settings can be specified to reflect patient-specific or experimental configurations.

Based on these inputs, OSS-DBSv2 computes the electric field in inhomogeneous and anisotropic brain tissue.
In subsequent post-processing steps, the software can estimate stimulation volumes and pathway activation for selected fiber tracts.

The sections below first provide a user-oriented guide and then detailed documentation of the main modeling components.
For more details about the first version of OSS-DBS, see [Butenko2019]_.

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   installation
   lead_dbs
   tutorial
   examples

.. toctree::
   :maxdepth: 1
   :caption: Modeling Components

   brain_geometry
   materials
   stimulation_signals
   electrodes
   point_analysis
   axon_models
   volume_conductor_model
   utils

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`

.. [Butenko2019] K. Butenko, C. Bahls, M. Schröder, R. Köhling and U. van Rienen,
   OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling,
   PLoS Comput Biol 16(7): e1008023.
   https://doi.org/10.1371/journal.pcbi.1008023
