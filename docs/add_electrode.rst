Adding an electrode
===================

Prepare the geometry
--------------------

* **Replace `electrodename`  by a descriptive name of the electrode that you want to add.**
* Create a file `electrodename_model.py` in `ossdbs/electrodes/electrode_models`.
* Implement a `dataclass` called `ElectrodenameParameters`, which shall contain all parameters
  required to build the electrode.
* Implement the actual electrode model in a class called `ElectrodenameModel`.
  A template for this class is given in :class:`ossdbs.electrodes.electrode_models.ElectrodeModel`.
* Create another file `electrodename.py` and fill it with a class called `Electrodename(ElectrodenameModel)`.
  This class does nothing else than building the electrode model with the default geometry.
  Customised electrode parameters are thus not possible using this class.

Final steps
-----------

* Import the electrode in the `__init__.py` file in the `ossdbs/electrodes/electrode_models` directory
* Add the electrode to `ossdbs/factories/custom_electrode_factory.py`
* Add the electrode to `ossdbs/factories/electrode_factory.py`
