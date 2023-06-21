# Test cases

Go to the respective directory and run then (except for input case 2)

```
ossdbs name_of_input_file.json
```

Output:

* Time signal at fixed points around the electrode
* Mesh
* Potential distribution at 130 Hz (`.vtk`)
* Impedance between active contact and ground


## Input case 1: homogeneous material

Boston Scientific Vercise electrode at X: 5, Y: 14, Z: -4.5.
The medium around is homogeneous.
One contact at 0.1 V and the other at 0 V.
EQS mode is active.

## Input case 2: UQ of SNEX 100 electrode parameters

Ten SNEX 100 geometries are simulated.
The electrode is placed at (X: 22.95, Y: 11.35, Z: 7.95) in the rat MRI `new_segmented_Atlas_GMWMCS`.
Contact 1 is at 0.2V and contact 2 at 0 V.
The encapsulation thickness is 0.1 mm and gray matter.
Simulation is run in FullSpectrum mode using QS at 130Hz and 60 us pulsewidth.

Runs with 

```
python3 electrode_uq.py input_case2.json

```

**TODO:**
Last test geometry currently fails.

## Input case 3: monopolar electrode

The monopolar electrode `MicroProbesCustomRodent` is implanted at (X: 22.95, Y: 11.35, Z: 7.95) in the rat MRI `new_segmented_Atlas_GMWMCS`.
Case grounding is used. **TODO: This contradicts the above point!** 
Simulation is run using QS at 130Hz and 60 us pulsewidth. **TODO: current-controlled mode with 200 muA not working!**
Results are calculated for one artificial neuron placed below the electrode tip.

## Input case 4: dielectric models

`AbbottStjudeActiveTip6146_6149` electrode is used at (X: 5, Y: 14, Z: -4.5).
The first contact is at 1 V and the second at 0 V.
The encapsulation layer is 0.3 mm thick and gray matter.
Simulation is run using QS at 130Hz and 60 us pulsewidth.
Instead of the conventional ColeCole4 model, we use constant values at 1 kHz for the permittivity and conductivity.
**TODO: Why permittivity when QS?**

## Input case 5


## Input case 6: Floating

**Do not run! Takes much too long and should be wrong.**

Standard human MRI (`icbm_avg_152_segmented.nii.gz`).
Position of the Boston Scientific Vercise elektroce: (X: 5, Y: 14, Z: -4.5).
Current-controlled stimulation with 1 mA on the lowest contact and ground at the contact above.
**TODO: What is the point of this configuration? Change!**

## Input case 7: Diffusion Tensor Image (DTI)

Use MRI and DTI data in MNI space from MNI_ICBM_2009b_NLIN_ASYM atlas from Lead-DBS templates.
Position of the Boston Scientific Vercise elektroce: (X: 5, Y: 14, Z: -4.5).
Voltage-controlled stimulation with 1V on the lowest contact and ground at the contact above.

**TODO: Implement processing of DTI data.**

## Input case 8: Point analysis for differenet areas

Use MRI data in MNI space from MNI_ICBM_2009b_NLIN_ASYM atlas from Lead-DBS templates.
Position of the Boston Scientific Vercise elektroce: (X: 5, Y: 14, Z: -4.5).
Voltage-controlled stimulation with 1V on the lowest contact and ground at the contact above.
We seed points in different patterns (line, array, or brain region) to evaluate the electric potential/field at different points.