Numpy style documentation
=========================
Use

.. code-block:: console

    $ sphinx-apidoc -f -o docs/ ../src

to generate new docstrings automaticly. Next step type

.. code-block:: console

    $ make html

in the docs folder.

Example Numpy Style documentation
---------------------------------
.. automodule:: src.example_numpy
    :members:
    :undoc-members:
    :show-inheritance:

