Adding an electrode
===================

Prepare the geometry
--------------------

* **Replace `electrodename`  by a descriptive name of the electrode that you want to add.**
* Create a file `electrodename_model.py` in `ossdbs/electrodes`.
* Implement a `dataclass` called `ElectrodenameParameters`, which shall contain all parameters
  required to build the electrode.
* Implement the actual electrode model in a class called `ElectrodenameModel`.
  A template for this class is given in :class:`ossdbs.electrodes.ElectrodeModel`.
  At the minimum, the following methods need to be implemented:
  `parameter_check`, `_construct_encapsulation_geometry`, `_construct_geometry`, `_body`, and `_contacts`.
* Import the model from the correct file to `ossdbs/electrodes/defaults.py` at the top of the file. 
  Also, add the default parameters of the electrode and a function named `Electrodename` to this file  
* Import the model from the correct file to `ossdbs/electrodes/__init__.py` and from `.defaults` at the top of the file.
  Additionally, add the model name/class name to the `ELECTRODE_MODELS`, `ELECTRODES`, `ELECTRODE_PARAMETERS`, and `__all__` dictionaries. 



