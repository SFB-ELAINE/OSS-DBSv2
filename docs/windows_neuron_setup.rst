NEURON Setup on Windows
========================

Follow these steps here to install ``NEURON`` on Windows.

Download NEURON
----------------

Download the official installer from:

- https://github.com/neuronsimulator/nrn/releases

Use the Windows ``.exe`` installer from the release you want.

Typical environment setup
--------------------------

NEURON needs to be reachable in two independent ways: the ``neuron`` Python
package must be on ``sys.path``, and Windows must be able to find the NEURON
DLLs (``libnrniv.dll`` etc) when ``neuron.hoc`` loads them.
This requires setting ``NEURONHOME``.

**Environment-wide** — simplest if you only use one Python environment:

.. code-block:: powershell

    $env:NEURONHOME = "C:\path\to\nrn"
    $env:PYTHONPATH = "C:\path\to\nrn\lib\python"
    $env:PATH = "C:\path\to\nrn\bin;$env:PATH"

Replace ``C:\path\to\nrn`` with the folder where NEURON was installed on your
machine (the installer's default is ``C:\nrn``).

**Scoped to one virtual environment** — preferred if you have multiple Python
installs/venvs on the machine, because a machine-wide ``PYTHONPATH``
impacts all virtual environments.
Instead, add it only to the venv that needs it via a ``.pth`` file,
and keep ``NEURONHOME`` set and
NEURON's ``bin`` directory on ``PATH`` persistently as a User
environment variable (``sysdm.cpl`` → Advanced → Environment Variables,
or in PowerShell):

.. code-block:: powershell

    [Environment]::SetEnvironmentVariable('NEURONHOME', 'C:\path\to\nrn', 'User')
    $userPath = [Environment]::GetEnvironmentVariable('PATH', 'User')
    [Environment]::SetEnvironmentVariable('PATH', "C:\path\to\nrn\bin;$userPath", 'User')

then set

.. code-block:: powershell

    "C:\path\to\nrn\lib\python" | Out-File -Encoding ascii .venv\Lib\site-packages\neuron.pth

``C:\path\to\nrn\bin`` must be on ``PATH``: OSS-DBS calls ``mknrndll``
to compile the NEURON model files for pathway activation modelling (PAM).
The ``import neuron`` check below works without it, so it will not catch
a missing ``PATH`` entry.

Verify that NEURON works
--------------------------

Run:

.. code-block:: powershell

    python -c "import neuron; print(neuron.__file__)"
    
or 

.. code-block:: powershell

    uv run python -c "import neuron; print(neuron.__file__)"

If this succeeds, the Python bindings are available in the current environment.
