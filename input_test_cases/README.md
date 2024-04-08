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

To conduct current-controlled stimulations, a fixed current for the stimulation amplitude can be provided in the settings. In this example, an electrode from Abbott St. Jude (`AbbottStJudeDirected6172`) is used with 0.1mA on the first and -0.1mA on the second contact. Also, current-controlled stimulation on multiple contacts is possible. To do so, we specify the amplitudes as follows: `C1: -3mA, C2: 1mA, C3: 1mA, C4: 1mA`.

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

WIP
