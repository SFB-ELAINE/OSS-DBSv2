.. OSS_DBS documentation master file, created by
   sphinx-quickstart on Fri Oct 21 11:27:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OSS-DBS's documentation!
===================================
.. note::
   This page is still under construction.

=========
Overview
=========
Deep Brain Stimulation} (DBS) is a widely used treatment for several motor and non-motor disorders.
Since the mechanism of action are not fully understood, computational models help us to predict the outcome of different treatments and optimize them.


OSS-DBS is a comprehensive tool to performe several DBS specific studies for humans, but also for animal studies in a highly automated workflow.
Therefore, the user can provide data like MRI and DTI data or selcet predefined datasets from the libary.
Also, the DBS electrodes can be selected from a predefined database which contains the most common electrode types, but could also be added manually.

The software performs calculations of the electric field within the inhomogenius and anisotropic brain tissue based on the given inputs.
Further, the activation of automaticly generated axons is clculated to estimate a volume of activated tissue (VAT) or specific realistic fiber tracts can be loaded to performe detailed pathway activation modeling (PAM).

A detailed overview over the implemented concepts can be find in the specific sections of this documentation.
For more details about the first veriosn of OSS-DBS you can also refer to [Butenko2019]_.

==============
Installation
==============
The software OSS-DBS can be easily installed using pip:

.. code-block:: console

    $ pip install oss_dbs

For first steps with OSS-DBS see the next section of this documentation.

=================
Table of contents
=================
.. toctree::
   :maxdepth: 1

   start
   examples
   volume_conductor_model
   electrodes


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. [Butenko2019] K. Butenko, C. Bahls, M. Schröder, R. Köhling and U. van Rienen, OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling, PLoS Comput Biol 16(7): e1008023. https://doi.org/10.1371/journal.pcbi.1008023.
