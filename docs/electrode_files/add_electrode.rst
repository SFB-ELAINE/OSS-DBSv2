Adding an Electrode
===================

Introduction
------------

This guide walks you through the process of adding a new electrode to OSS-DBS. 
Follow these steps to ensure seamless integration of your custom electrode.

Prepare the Geometry
--------------------

**Before you begin, replace `electrodename` with a descriptive name for your electrode in the following steps.**

.. note:: When following the steps below, you can use an existing electrode as a template and modify it to suit your new design.

1. **Create a New File:**  
   Create a file named `electrodename_model.py` within the `ossdbs/electrodes` directory.

2. **Define Electrode Parameters:**  
   Implement a `dataclass` named `ElectrodenameParameters`. This class should include all parameters required for constructing the electrode geometry and defining its properties.

3. **Implement the Electrode Model:**  
   Create a class named `ElectrodenameModel` to define the electrode model. Use the base class :class:`ossdbs.electrodes.ElectrodeModel` as a template. At a minimum, ensure the following methods are implemented:
   - `parameter_check`  
   - `_construct_encapsulation_geometry`  
   - `_construct_geometry`  
   - `_body`  
   - `_contacts`  
   - `get_center_first_contact`  
   - `get_distance_l1_l4`  

4. **Add Default Parameters:**  
   Import the model into the `ossdbs/electrodes/defaults.py` file. At the beginning of this file:
   - Define the default parameters for your electrode in a `dataclass` or dictionary.  
   - Add a function named `Electrodename()` to return the electrode's default parameter instance.  

5. **Update Initialization File:**  
   Import the new module into the `ossdbs/electrodes/__init__.py` file. At the top of the file:
   - Add the module to the imports list.  
   - Include your model in the following dictionaries:
   
     - `ELECTRODE_MODELS`  
     - `ELECTRODES`  
     - `ELECTRODE_PARAMETERS`  
     - `__all__`  

Following these steps ensures that the new electrode is fully integrated into the OSS-DBS framework and is available for use.

Tips for Testing
-----------------

- Use the provided templates and methods to verify the geometry and parameters of your electrode.
- Ensure that the parameter names, data types, and default values are consistent with the conventions used for existing electrodes.
- Implement unit tests to validate the new electrode's behavior and ensure compatibility with the rest of the software.

Final Verification
-------------------

Once the new electrode is integrated:
1. Generate the Sphinx documentation for your electrode and verify its appearance.
2. Check that the electrode appears correctly in the list of available electrodes in OSS-DBS.
3. Validate its geometry and functionality through simulation or testing with existing workflows.
