Human DBS Example
=================

This example demonstrates DBS modeling on a human subject using MRI- and DTI-based tissue characterization and a Boston Scientific Vercise directed electrode.  
It illustrates how to configure anatomical data, electrode orientation and placement, dielectric modeling, solver settings, and typical output options.

This example is representative of patient-specific DBS simulations based on real neuroimaging data.


Input Highlights
----------------

- **BrainRegion:** spherical subvolume around the stimulation site, sufficiently large to include local tissue heterogeneity and reduce boundary effects.
- **Electrode:** Boston Scientific Vercise directed lead with its true shaft orientation in the patient's anatomical coordinate system.
- **MRI:** used to distinguish CSF, GM, WM, and blood.
- **DTI:** provides local anisotropy (tensor conductivity).
- **Dielectric model:** 4-Cole-Cole dispersion model (standard for brain tissue).
- **Mesh:** adaptive strategy with fine resolution near the stimulation region.
- **Stimulation:** frequency-domain solution at 10 kHz using a multisine formulation.
- **Outputs:** impedance estimate, electrode geometry, and field/VTP data exported for visualization and post-analysis.


Full Input Example
------------------

.. code-block:: json

    {
      "BrainRegion": {
        "Center": {"x[mm]": -13.99, "y[mm]": -7.73, "z[mm]": -7.91},
        "Dimension": {"x[mm]": 70.0, "y[mm]": 70.0, "z[mm]": 70.0},
        "Shape": "Sphere"
      },
      "Electrodes": [
        {
          "Name": "BostonScientificVerciseDirected",
          "CustomParameters": null,
          "Rotation[Degrees]": 0.0,
          "Direction": {"x[mm]": -0.45, "y[mm]": 0.65, "z[mm]": 0.61},
          "TipPosition": {"x[mm]": -13.99, "y[mm]": -7.73, "z[mm]": -7.91}
        }
      ],
      "MaterialDistribution": {
        "MRIPath": "./input_files/sub-John_Doe/JD_segmask.nii.gz",
        "MRIMapping": {
          "Unknown": 0,
          "CSF": 3,
          "White matter": 2,
          "Gray matter": 1,
          "Blood": 4
        },
        "DiffusionTensorActive": true,
        "DTIPath": "./input_files/sub-John_Doe/JD_DTI_NormMapping.nii.gz"
      },
      "DielectricModel": {"Type": "ColeCole4", "CustomParameters": null},
      "Mesh": {
        "LoadMesh": false,
        "LoadPath": "",
        "MeshingHypothesis": {
          "Type": "Fine",
          "MaxMeshSize": 1000.0,
          "MeshSizeFilename": ""
        },
        "SaveMesh": false
      },
      "StimulationSignal": {
        "CurrentControlled": false,
        "Type": "Multisine",
        "ListOfFrequencies": [10000.0]
      },
      "Solver": {
        "Type": "CG",
        "Preconditioner": "bddc",
        "PreconditionerKwargs": {},
        "PrintRates": true,
        "MaximumSteps": 200,
        "Precision": 1e-8
      },
      "PointModel": {
        "Pathway": {
          "Active": false,
          "FileName": ""
        },
        "Lattice": {
          "Active": true,
          "Center": {"x[mm]": -13.99, "y[mm]": -7.73, "z[mm]": -7.91},
          "Shape": {"x": 20, "y": 20, "z": 20},
          "Direction": {"x[mm]": -0.45, "y[mm]": 0.65, "z[mm]": 0.61},
          "PointDistance[mm]": 0.5
        }
      },
      "OutputPath": "Results",
      "ComputeImpedance": true,
      "ExportVTK": true,
      "ExportElectrode": true
    }


Explanation of Key Parameters
-----------------------------

**BrainRegion:**  
A 70 × 70 × 70 mm spherical region centered at the targeted stimulation site.
This ensures:
- inclusion of relevant neuroanatomy,
- realistic current spread,
- minimized artificial boundary effects.

**Electrode Orientation:**  
The ``Direction`` vector aligns the shaft with the patient's anatomy.
This supports:
- directional stimulation,
- contact-wise modeling,
- and consistent mapping with tractography.

**Tissue Properties (MRI + DTI):**
- MRI segmentation defines local tissue class assignment (GM, WM, CSF…).
- DTI maps anisotropic conductivity tensors into OSS-DBSv2.

**Dielectric Model:**  
Cole-Cole dispersion (4-term) is the standard frequency-dependent model for brain tissue.  
No custom override is used here to preserve realistic electrical behavior.

**Meshing:**  
``Fine`` hypothesis ensures adequate resolution in the vicinity of the electrode.
The mesh is generated automatically.

**Stimulation:**  
- Frequency-domain solution at 10 kHz.
- Multisine allows solving a single frequency slice efficiently.

**Outputs:**  
- Impedance provides model diagnostics and electrode–tissue interface estimation.
- VTK files allow visualization (ParaView, Lead-DBS, MNE, etc.).
- The electrode geometry is exported for debugging and high-quality plotting.


Running the Example
-------------------

Execute the simulation from the command line:

.. code-block:: bash

    ossdbs input_human.json

A logfile and the result files (VTK, impedance, and optional outputs) will be written to the ``Results/`` folder.


Typical Uses
------------

This type of patient-specific setup is used for:

- evaluating contact selection,
- comparing different stimulation parameters,
- impedance analysis,
- grounding and reference validation,
- tract-related field interpretation,
- pre/post clinical electrode position review.
