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
  The number of contacts and two methods need to be implemented:
  `_construct_geometry` and `_construct_encapsulation_geometry`.
* Add the default parameters of the electrode to `ossdbs/electrodes/defaults.py`. 
  Add a function to `ossdbs/electrodes/__init__.py` named `Electrodename`.
  This function does nothing else than building the electrode model with the default geometry.
  Customised electrode parameters are thus not possible using this class.

Final steps
-----------

* Add the electrode to `ossdbs/factories/custom_electrode_factory.py`
* Add the electrode to `ossdbs/factories/electrode_factory.py`
