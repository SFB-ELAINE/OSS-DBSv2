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
Deep Brain Stimulation (DBS) is a widely used treatment for several motor and non-motor disorders.
Since the mechanism of action are not fully understood, computational models help us to predict the
outcome of different treatments and optimize them.


OSS-DBS is a comprehensive tool to performe several DBS specific studies for humans, but also for animal studies in a highly automated workflow.
Therefore, the user can provide data like MRI and DTI data or use the data from the example files.
Further, there is a set of predifined DBS electrodes which contains the most common electrode types,
but you can add addidtional electrodes easily.

The software performs calculations of the electric field within the inhomogenius and anisotropic brain
tissue based on the given data and parameters.

In the steps of postprocessing the software calculates the volume of activated tissue (VAT) or the
activation of specific axons and realistic fiber tracts to performe detailed pathway activation
modeling (PAM).

A detailed overview over the implemented concepts can be find in the related sections of this
documentation. For more details about the first version of OSS-DBS you can refer to [Butenko2019]_.

============
Installation
============
In near future, the software OSS-DBS can be easily installed using pip:

.. code-block:: console

    $ pip install ossdbs

For now, please use:

.. code-block:: console

    $ git clone https://github.com/SFB-ELAINE/OSS-DBSv2.git 
    $ pip install . 


For first steps with OSS-DBS see the next section of this documentation.
Some examples to run you can find in the :ref:`Examples <examples>` section.

For Developers
==============

Install 

.. code-block:: console

   $ pip install -e ".[dev]"

Or if you want to update the Sphinx documentation locally install

.. code-block:: console

   $ pip install -e ".[all]"

=================
Table of contents
=================
.. toctree::
   :maxdepth: 1

   start
   examples
   brain_geometry
   materials
   stimulation_signals
   electrodes
   point_analysis
   axon_models
   volume_conductor_model
   utils


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. [Butenko2019] K. Butenko, C. Bahls, M. Schröder, R. Köhling and U. van Rienen, OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling, PLoS Comput Biol 16(7): e1008023. https://doi.org/10.1371/journal.pcbi.1008023.
