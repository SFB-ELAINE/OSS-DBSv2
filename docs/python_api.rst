Python scripting API
====================

OSS-DBSv2 can be driven entirely from Python without the CLI. This is useful
for parameter sweeps, custom post-processing, embedding OSS-DBSv2 in larger
pipelines, or interactive exploration in notebooks.

Quick start
-----------

The simplest way to run a simulation from Python is to load a JSON settings
dictionary and call the same pipeline that the CLI uses:

.. code-block:: python

   import json
   import ossdbs

   ossdbs.set_logger()

   with open("input.json") as f:
       settings = json.load(f)

   from ossdbs.main import main_run
   main_run(settings)

This is equivalent to ``ossdbs input.json`` on the command line.

Step-by-step API
----------------

For more control, individual pipeline steps can be called separately. The
typical sequence mirrors the CLI pipeline:

.. code-block:: python

   import json
   import ngsolve
   import ossdbs
   from ossdbs.api import (
       build_brain_model,
       generate_electrodes,
       load_images,
       prepare_dielectric_properties,
       prepare_solver,
       prepare_stimulation_signal,
       prepare_volume_conductor_model,
       run_volume_conductor_model,
       set_contact_and_encapsulation_layer_properties,
       validate_solver_settings,
   )
   from ossdbs.fem import ConductivityCF
   from ossdbs.model_geometry import ModelGeometry
   from ossdbs.utils.settings import Settings

   ossdbs.set_logger()

   # 1. Load and complete settings
   with open("input.json") as f:
       input_settings = json.load(f)
   settings = Settings(input_settings).complete_settings()

   # 2. Load imaging data
   mri_image, dti_image = load_images(settings)

   # 3. Build geometry
   electrodes = generate_electrodes(settings)
   brain_model = build_brain_model(settings, mri_image)
   geometry = ModelGeometry(brain_model, electrodes)
   set_contact_and_encapsulation_layer_properties(settings, geometry)

   # 4. Prepare material model
   validate_solver_settings(settings, geometry)
   dielectric_properties = prepare_dielectric_properties(settings)
   materials = settings["MaterialDistribution"]["MRIMapping"]
   conductivity = ConductivityCF(
       mri_image,
       brain_model.brain_region,
       dielectric_properties,
       materials,
       geometry.encapsulation_layers,
       complex_data=settings["EQSMode"],
       dti_image=dti_image,
   )

   # 5. Create solver and volume conductor
   with ngsolve.TaskManager():
       solver = prepare_solver(settings)
       volume_conductor = prepare_volume_conductor_model(
           settings, geometry, conductivity, solver
       )
       frequency_domain_signal = prepare_stimulation_signal(settings)

       # 6. Refine mesh and solve
       volume_conductor.prepare_mesh_refinements(
           settings["Mesh"]["MaterialRefinementSteps"]
       )
       run_volume_conductor_model(
           settings, volume_conductor, frequency_domain_signal
       )

After the solve, the volume conductor object provides access to the solution:

.. code-block:: python

   # Access the potential at the last solved frequency
   potential = volume_conductor.potential

   # Evaluate at arbitrary points (Nx3 array in mm)
   import numpy as np
   points = np.array([[0.0, 0.0, 5.0], [0.0, 0.0, 6.0]])
   values = volume_conductor.evaluate_potential_at_points(points)

   # Impedance data (if ComputeImpedance was enabled)
   impedances = volume_conductor.impedances

Key functions
-------------

The ``ossdbs.api`` module provides the building blocks:

- ``load_images(settings)`` — load MRI and optional DTI data
- ``generate_electrodes(settings)`` — create electrode CAD models
- ``build_brain_model(settings, mri_image)`` — create the brain geometry
- ``set_contact_and_encapsulation_layer_properties(settings, geometry)`` —
  apply contact activation, voltages, currents, and surface impedance
- ``prepare_dielectric_properties(settings)`` — set up the dielectric model
- ``prepare_solver(settings)`` — instantiate solver and preconditioner
- ``prepare_volume_conductor_model(settings, geometry, conductivity, solver)``
  — create the appropriate volume conductor variant
- ``prepare_stimulation_signal(settings)`` — build the frequency-domain signal
- ``run_volume_conductor_model(settings, vc, signal)`` — solve at all
  frequencies, export results
- ``generate_point_models(settings)`` — create lattice or pathway evaluators

Settings handling
-----------------

The ``Settings`` class in ``ossdbs.utils.settings`` merges user-provided
values with defaults. Call ``Settings(input_dict).complete_settings()`` to
get a fully populated dictionary before passing it to any API function.
This is the same step the CLI performs internally.

Related pages
-------------

- :doc:`tutorial` — first simulation from the command line
- :doc:`input_settings` — JSON input reference (all keys and defaults)
- :doc:`volume_conductor_model` — solver, mesh, and formulation details
- :doc:`stimulation_signals` — signal types and frequency-domain settings
- :doc:`examples` — runnable example configurations

API reference
-------------

.. automodule:: ossdbs.api
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.main
    :members:
    :undoc-members:
    :show-inheritance:
