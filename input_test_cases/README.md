OSS-DBS v2.0 Test Cases
=======================

This folder contains various test cases showcasing the functionality of OSS-DBS v2.0. To execute different cases, use the command line as follows:

```bash
ossdbs path_to_folder/input_file.json
```

The command line displays a dialog while the software is running. Upon completion, all outputs are stored in a new folder, which is, by default, inside the folder where the input.json is stored.

The output folder contains different files, depending on the chosen studies. Also, the exporting of various files can be determined at the end of the input.json file. The following results can be stored:

### Results:

* Material and conductivity distribution as `.vtk` file,
* Electric potential and field distribution as `.vtk` file,
* Probed electric potential and field as `.h5` file,
* Volume of tissue activated as `.nii` file,
* Used electrodes as `.vtk` file,
* Create mesh as `.vol.gz` file,
* Estimated impedance as `.csv` file.

Case 1: Brain Material
----------------------

Our first test case contains two input dictionaries to demonstrate the use of inhomogeneous and anisotropic tissue properties. The first one uses a homogeneous MRI image but no diffusion tensor image (DTI), whereas the other input file uses an MRI scan of a human brain.
A Boston Scientific Vercise electrode (`BostonScientificVercise`) is used and implanted at the subthalamic nucleus (STN). A unit amplitude with 1V on the lowest contact and grounding on the second contact is modeled. The tissue properties are estimated based on the ColeCole4 model at a single frequency of 10kHz. The quasi-static (QS) approximation of Maxwell's equation is solved, and the outputs are stored in the results folder.

Case 2: Custom Parameters
-------------------------

The second test case shows how to use custom parameters for electrode geometries or material models. Here, the directed Boston Scientific Vercise electrode (`BostonScientificVerciseDirected`) is placed at the STN with an encapsulation layer of 0.2 mm surrounding the electrode. The encapsulation layer is assumed to consist of gray matter but can be changed to any provided tissue type. In the `input_custom_electrode.json`, we demonstrate how to modify an electrode by slightly varying the length of the contacts. In the second input dictionary, we change the material model to a constant value for all tissue types, where the corresponding values for the conductivity and permittivity are defined in the inputs.

Case 3: Case Grounding
----------------------

In this case, we use a custom-designed monopolar electrode by MicroProbes (`MicroProbesRodentElectrode`), which is used for DBS in rodents. The electrode is implanted in the STN within the same MRI as used before. The stimulation amplitude of the only contact is 1V, and the outer boundary of the brain region is assumed to be the grounding. In the case of a custom brain shape, specified regions can also be used as grounding in the same way.


Case 4: Current Controlled Stimulation
--------------------------------------

To conduct current-controlled stimulations, a fixed current for the stimulation amplitude can be provided in the settings. In `input_current_controlled.json`, an electrode from Abbott St. Jude (`AbbottStJudeDirected6172`) is used with 0.1mA on the first and -0.1mA on the second contact. Current-controlled stimulation on multiple contacts is demonstrated in `input_multi_current.json`, where four contacts are driven simultaneously with `C1: -3mA, C2: 1mA, C3: 1mA, C4: 1mA`. The per-contact (floating) potentials resulting from such a multi-contact configuration are written to `floating_potentials.csv`.

Case 5: Stimulation Signals
---------------------------

The Medtronic SenSight electrode (`MedtronicSenSightB33005`) is used to simulate a stimulation with a rectangular signal with a base frequency of 130 Hz and 60 us pulse width. After transforming the time signal into the frequency domain and solving in the frequency domain, the results are transferred back into the time domain. In this test case, the OcatveBand method is used to reduce the number of computed frequencies.

Case 6: Floating Contacts
-------------------------

To demonstrate the use of floating contacts, a new `.json` is created using the PINS Medical electrode (`PINSMedicalL303`). The first contact uses a 1V stimulation amplitude, and the third contact is used as ground. The second contact is modeled as a floating conductor.

Case 7: Volume of Tissue Activated (VTA)
----------------------------------------

The contained input dictionary uses a uniform grid around the active contact of the DIXI SEEG electrode (`DixiSEEG10`) to estimate the electric field at those points and threshold it by a specified value. As a result, the estimated VTA is stored in Nifty format. As an additional option, the electrode can virtually be removed from the VTA to collapse it inside the electrode's trajectory.

Case 8: Pathway Activation Modelling (PAM)
------------------------------------------

Using a DBS electrode from Medtronic (`Medtronic3387`) placed in the STN, the electric field is evaluated along the points of an axon. Therefore, the coordinates are given in a structured `.h5`  file, and the path to the file is provided in the inputs. The `.h5` file needs to be structured in groups (e.g., pathways), where each group can contain multiple datasets (e.g., axons). Each axon includes an array of 3D coordinates and needs to follow the naming convention (axon0, axon1, axon2, etc.). The results are stored in the output folder in the same structure in `.h5` format.

Case 9: StimSets
----------------

Run pathway activation modeling (PAM) for multiple stimulation protocols in a single FEM solve. The `Medtronic3387` electrode is used together with the same axon `.h5` file as in Case 8. Each contact is solved once in unit-current mode, and the resulting per-contact fields are then linearly combined according to the rows of a stimulation-set table provided as `Current_protocol.csv`. Each row of that CSV lists the current (in mA) applied to every contact for one protocol; the PAM step is then repeated for every row. Results for each contact's unit solution are written to `Results_PAME1C1`, `Results_PAME1C2`, … (one folder per contact), while the pathway-activation outcomes for each protocol are written alongside the main `Results_PAM` folder.

Case 10: Floating Contacts with Surface Impedance
-------------------------------------------------

This case extends Case 6 to contacts with a lumped surface impedance (Robin boundary condition), using the `PINSMedicalL303` electrode. The following variants are provided:

* `input_floating.json`: baseline, identical in spirit to Case 6 (1V active contact, floating middle contact, grounded third contact, no surface impedance).
* `input_floating_encap.json`: same configuration but with a 0.2 mm encapsulation layer around the electrode.
* `input_floating_impedance.json`: the floating middle contact is assigned a surface impedance (resistive model `R = 1e3 Ohm*mm^2`) to mimic a non-ideal electrode-tissue interface. `HPRefinement` is enabled.
* `input_floating_impedance_noground.json`: all three contacts are floating and current-controlled (C1: +1 mA, C3: -1 mA, C2: 0 mA), each with its own surface impedance. The computed per-contact floating potentials and currents are written to `floating_potentials.csv` / `currents.csv`.

The helper script `plot_comparison.py` produces a side-by-side comparison of the different variants.

Case 11: Surface Impedance (Lempka2009 models)
----------------------------------------------

This case demonstrates the difference between interface impedance modeled as a lumped equivalent circuit and surface impedance modeled as a per-contact Robin boundary condition. A `BostonScientificVercise` electrode is placed in a homogeneous medium, and the stimulation signal is a multisine sweep from 1 Hz to 10 kHz so that the impedance spectrum can be evaluated. The provided input files are:

* `input_no_interface.json`: no electrode–tissue interface (baseline).
* `input_homogeneous.json` / `input_homogeneous_lempka2009.json`: interface impedance via a lumped circuit.
* `input_homogeneous_surface_impedance_lempka2009.json` / `input_surface_impedance_lempka2009.json`: per-contact surface impedance using the CPE double-layer model of Lempka et al. (2009), i.e. `"SurfaceImpedance": {"Model": "CPE_dl", "Parameters": {"dl_k": 1.5e6, "dl_alpha": 0.8}}`.

The helper script `get_mesh.py` can be used to generate the shared mesh, and `_impedance_overview.pdf` shows the resulting impedance spectra for the different models.

Running the test suite
----------------------

The cases above are wired into a pytest-based runner (`test_simulations.py`). Tests are marked with `simulation`, and with additional marks (`slow`, `requires_neuron`, `vta`, `pam`, `surface_impedance`) so that subsets can be selected:

```bash
# run everything
pytest input_test_cases/test_simulations.py -m simulation -v

# fast cases only
pytest input_test_cases/test_simulations.py -m "simulation and not slow" -v

# a specific case by its ID
pytest input_test_cases/test_simulations.py -k "brain_material" -v
```

Each test runs the simulation (and, where applicable, the PAM step) and compares the generated outputs against the baselines stored under `desired_output/`. Baselines can be refreshed with `update_desired_output.py`.
