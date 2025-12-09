Installation
============

OSS-DBS can now be easily installed using pip.
The software runs with Python versions 3.10, 3.11, and 3.12.

Basic Installation
------------------

For most users, simply run the following command:

.. code-block:: console

    $ pip install ossdbs

On macOS and Windows, the NEURON simulator must be installed separately **before** installing OSS-DBSv2.
Installation instructions for NEURON are available here: [https://nrn.readthedocs.io/en/latest/install/install.html].

Developer Installation
----------------------

To install OSS-DBSv2 from source for development:

.. code-block:: console

    $ git clone https://github.com/SFB-ELAINE/OSS-DBSv2.git
    $ cd OSS-DBSv2
    $ pip install -e .

Additional Installation Options
--------------------------------

- Install OSS-DBSv2 together with the test suite:

  .. code-block:: console

      $ pip install -e ".[test]"

- Install OSS-DBSv2 with development dependencies:

  .. code-block:: console

      $ pip install -e ".[dev]"

- Install OSS-DBSv2 with *all* optional dependencies:

  .. code-block:: console

      $ pip install -e ".[all]"


Next Steps
----------

For first steps with OSS-DBSv2, see the :ref:`Tutorial <tutorial>` section.
Example simulations can be found in :ref:`Examples <examples>`.