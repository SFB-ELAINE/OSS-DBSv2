Installation
============

OSS-DBS can now be easily installed using pip.
The software runs with Python versions 3.11, 3.12, 3.13, and 3.14.

Basic Installation
------------------

For most users, simply run the following command:

.. code-block:: console

    $ pip install ossdbs

On Windows, the NEURON simulator must be installed separately **before** installing OSS-DBSv2,
since NEURON has no PyPI wheel for Windows. See :doc:`windows_neuron_setup` for
step-by-step instructions, including how to make NEURON importable from a
virtual environment and how to set ``NEURONHOME``.

Developer Installation
----------------------

To install OSS-DBSv2 from source for development:

.. code-block:: console

    $ git clone https://github.com/SFB-ELAINE/OSS-DBSv2.git
    $ cd OSS-DBSv2
    $ pip install -e .

Additional Installation Options
-------------------------------

- Install OSS-DBSv2 together with the test suite:

  .. code-block:: console

      $ pip install -e ".[test]"

- Install OSS-DBSv2 with development dependencies:

  .. code-block:: console

      $ pip install -e ".[dev]"

- Install OSS-DBSv2 with *all* optional dependencies:

  .. code-block:: console

      $ pip install -e ".[all]"


Checking the installation
-------------------------

After installation, a simple first check is to run one of the shipped example
input files as described in :doc:`tutorial`.

Troubleshooting: ``GLIBCXX_3.4.32`` not found (Linux)
-----------------------------------------------------

Since NEURON 9, the axon models (MOD files) are compiled to C++ with the
compiler installed on your system the first time a pathway activation is
computed. A library compiled with GCC 13 or newer needs a C++ runtime
(libstdc++) that provides ``GLIBCXX_3.4.32``. If your Python environment
brings an older libstdc++, typically a conda environment such as the one Lead-DBS
creates, loading the compiled models fails with an error like:

.. code-block:: console

    version `GLIBCXX_3.4.32' not found (required by .../neuron_model/x86_64/libnrnmech.so)

Update the C++ runtime of the environment so that it is at least as new as your
system compiler, for example:

.. code-block:: console

    $ conda activate <your-ossdbs-environment>
    $ conda install -c conda-forge "libstdcxx>=13"

If you use a compiler that is not the system default (e.g., from an HPC module),
make sure its libstdc++ is also found when you run OSS-DBSv2, for example by
loading the same module or adding its library folder to ``LD_LIBRARY_PATH``.

Next Steps
----------

For first steps with OSS-DBSv2, see :doc:`tutorial`.
Example simulations can be found in :doc:`examples`.
