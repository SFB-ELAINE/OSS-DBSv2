Installation
============

OSS-DBS can now be easily installed using pip. The software is tested with Python versions 3.8, 3.9, 3.10, and 3.11.

Basic Installation
------------------

For most users, simply run the following command:

.. code-block:: console

    $ pip install ossdbs

Windows users and Mac users with Python 3.8: Please install NEURON separately before installing OSS-DBS. The instructions can be found here: [https://nrn.readthedocs.io/en/latest/install/install.html].

Developer Installation
----------------------

To install OSS-DBS for development purposes, clone the repository into a local directory, navigate to the directory, and run:

.. code-block:: console

    $ git clone https://github.com/SFB-ELAINE/OSS-DBSv2.git
    $ cd OSS-DBSv2
    $ pip install -e .

Additional Installation Options
--------------------------------

- To run the test suite:

  .. code-block:: console

      $ pip install -e ".[test]"

- To develop OSS-DBS:

  .. code-block:: console

      $ pip install -e ".[dev]"

- To install all of the above:

  .. code-block:: console

      $ pip install -e ".[all]"

For first steps with OSS-DBS, see the next section of this documentation. Some examples to run can be found in the :ref:`Examples <examples>` section.