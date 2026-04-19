.. _electrodes:

Electrodes
==========

This is an overview about the electrode models used for simulation. Each
electrode model defines the geometry (tip, contacts, shaft) and default
dimensions. The electrode is placed in the simulation domain via
``TipPosition``, ``Direction``, and ``Rotation[Degrees]`` in the JSON input.

Custom electrode parameters
---------------------------

Every electrode model has a corresponding parameter class that defines its
geometric dimensions (tip length, contact length, contact spacing, lead
diameter, total length, etc.). To modify these dimensions, append ``Custom``
to the electrode name and provide a ``CustomParameters`` dictionary:

.. code-block:: json

   {
     "Name": "BostonScientificVerciseDirectedCustom",
     "CustomParameters": {
       "tip_length": 1.5,
       "contact_length": 1.5,
       "contact_spacing": 0.5,
       "lead_diameter": 1.3,
       "total_length": 450.0
     },
     "Rotation[Degrees]": 0.0,
     "Direction": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0},
     "TipPosition": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0}
   }

The available parameter names depend on the electrode model. Common parameters
across most models include ``tip_length``, ``contact_length``,
``contact_spacing``, ``lead_diameter``, and ``total_length`` (all in mm).
See the API reference for each electrode's parameter class.

.. note::

   When using the standard (non-Custom) electrode name, default dimensions are
   used and ``CustomParameters`` is ignored.

Available electrode models
--------------------------

.. toctree::
   :maxdepth: 1

   electrode_files/Abbott_Active_Tip
   electrode_files/Abbott_StJude_Directed
   electrode_files/Boston_Scientific_Vercise_Directed
   electrode_files/Boston_Scientific_Vercise
   electrode_files/Boston_Scientific_Cartesia
   electrode_files/Dixi_Microtechniques
   electrode_files/Medtronic_DBS
   electrode_files/Medtronic_SenSight
   electrode_files/MicroElectrode
   electrode_files/Micro_Probes_Custom_Rodent
   electrode_files/Micro_Probes_Snex100
   electrode_files/Neuro_Pace
   electrode_files/Pins_Medical
   electrode_files/PMTsEEG
   electrode_files/add_electrode

.. note:: The total length does not influence the computational domain; therefore, the above electrodes are modeled at 400 mm.

Related pages
-------------

- :doc:`input_settings` — ``Electrodes`` JSON configuration and contact settings
- :doc:`brain_geometry` — how electrodes are integrated into the model geometry
- :doc:`volume_conductor_model` — boundary conditions on contacts, including
  :ref:`surface impedance <surface-impedance>` and
  :ref:`stimulation modes <stimulation-modes>`