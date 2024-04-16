Adding an electrode
===================

Introduction
------------

This guide walks you through the process of adding a new electrode to OSS-DBS. 
Follow these steps to ensure seamless integration of your custom electrode.

Prepare the geometry
--------------------

**Before you begin, remember to replace electrodename with a descriptive name for your electrode in the following steps.**

.. note::  When following the steps below, you can copy an existing electrode and modify it according to your new design.

1. Create a file named `electrodename_model.py` within the `ossdbs/electrodes` directory.
2. Implement a `dataclass` named `ElectrodenameParameters`, which encompass all parameters necessary for constructing the electrode.
3. Implement the electrode model within a class named `ElectrodenameModel`. You can utilize the template provided in :class:`ossdbs.electrodes.ElectrodeModel`. Ensure that the following methods are implemented at a minimum: `parameter_check`, `_construct_encapsulation_geometry`, `_construct_geometry`, `_body`, `_contacts`, `get_center_first_contact`, and `get_distance_l1_l4`.
4. Import the model into the `ossdbs/electrodes/defaults.py` file at the beginning of the file. Additionally, add the default parameters of the electrode and a function named `Electrodename` to this file. 
5. Import the module into the `ossdbs/electrodes/__init__.py` file, along with importing from the `.defaults` at the top of the file. Furthermore, add the model name/class name to the `ELECTRODE_MODELS`, `ELECTRODES`, `ELECTRODE_PARAMETERS`, and `__all__` dictionaries. 

Following these steps ensures proper integration of the new electrode into the software environment.
