Benchmarking OSS-DBS v.0.1
===========================

.. note::
    Need to find a good marker for benchmarking!

This Benchmark with OSS-DBS v.0.1 [Butenko2019]_ will show that results from the new version are consistent
with previous results. Also we will the improvements in the handling of the new version as well as the advanteges
of the performance.

Benchmarking case
------------------

For simplification we show the benchmark for an inhomogenius but isotropic case by using the quasi-static (QS)
approximation of Maxwells' equation. Therefore we use a segmented ICBM 152 human brain atlas [Fonov2011]_ and place an electrode in
the human subthamalic nucleus (STN). The electrode is specified as Medtronic 3389 with a voltage of 1V at the
lowest contact and the third contact is used as ground contact.

Input
------

Here you can see how the input file in the new verion:


.. code-block:: bash

     "Electrodes": [
        {
            "Name": "MicroProbesCustomRodent",
            "Rotation": 0.0,
            "Direction": [0.0, 0.0, 1.0],
            "Translation": [5, 5, 5],
            "Contacts": {
                "Active": [true],
                "Value": [1.0]
            }
        }
    ],
    "FEMMode": "QS",
    "MagneticResonanceImage":
        {
            "MaterialCoding": {
                "Unknown": 0,
                "GrayMatter": 1,
                "WhiteMatter": 2,
                "CerebrospinalFluid": 3
            },
            "Path": "./input_files/TestMRI.nii"
        },
    "MeshElementOrder": 2,
    "OutputPath": "test_result",
    "RegionOfInterest":
        {
            "Active": false,
            "Center": [0, 0, 0],
            "Shape": [1, 1, 1]
        },
    "SpectrumMode": "NoTruncation",
    "StimulationSignal":
        {
            "Type": "Rectangle",
            "Frequency": 130.0,
            "PulseWidthPercentage": 0.0078,
            "TopWidthPercentage": 0.0
        }


The coresponding input file to the old version has the following form:

.. code-block:: bash

    d = {
    'voxel_arr_MRI': 0,
    'voxel_arr_DTI': 0,
    'Init_neuron_model_ready': 0,
    'Init_mesh_ready': 0,
    'Adjusted_neuron_model_ready': 0,
    'CSF_mesh_ready': 0,
    'Adapted_mesh_ready': 0,
    'signal_generation_ready': 0,
    'Parallel_comp_ready': 0,
    'Parallel_comp_interrupted': 0,
    'IFFT_ready': 0,
    'MRI_data_name': 'icbm_avg_152_segmented.nii.gz',
    'MRI_in_m': 0,
    'DTI_data_name': '',
    'DTI_in_m': 0,
    'CSF_index': 1.0,
    'WM_index': 3.0,
    'GM_index': 2.0,
    'default_material': 3,
    'Electrode_type': 'Medtronic3389',
    'Brain_shape_name': '0',
    'x_length': 30.0,
    'y_length': 30.0,
    'z_length': 30.0,
    'Aprox_geometry_center': [10.92957028, -12.11697637, -7.69744601],
    'Implantation_coordinate_X': 10.929,
    'Implantation_coordinate_Y': -12.117,
    'Implantation_coordinate_Z': -7.697,
    'Second_coordinate_X': 10.929,
    'Second_coordinate_Y': -9.437,
    'Second_coordinate_Z': 3.697,
    'Rotation_Z': 0.0,
    'encap_thickness': 0.3,
    'encap_tissue_type': 2,
    'encap_scaling_cond': 0.5,
    'encap_scaling_perm': 1.0,
    'pattern_model_name': '0',
    'diam_fib': [5.7],
    'n_Ranvier': [35],
    'v_init': -80.0,
    'Neuron_model_array_prepared': 0,
    'Name_prepared_neuron_array': 'Vert_for_VTA_SRI_space.csv',
    'Global_rot': 1,
    'x_seed': 10.929,
    'y_seed': -12.117,
    'z_seed': -7.697,
    'x_steps': 9,
    'y_steps': 0,
    'z_steps': 9,
    'x_step': 1.0,
    'y_step': 1.0,
    'z_step': 1.0,
    'alpha_array_glob': [0],
    'beta_array_glob': [0],
    'gamma_array_glob': [0],
    'X_coord_old': 0,
    'Y_coord_old': 0,
    'Z_coord_old': 0,
    'YZ_angles': [0],
    'ZX_angles': [0],
    'XY_angles': [0],
    'EQS_core': 'EQS',
    'Skip_mesh_refinement': 1,
    'refinement_frequency': [0],
    'num_ref_freqs': -1,
    'rel_div_CSF': -1,
    'Adaptive_frac_div': 0.0,
    'Min_Scaling': 1.0,
    'CSF_ref_reg': 0.0,
    'rel_div': 0.0,
    'rel_div_current': 0.0,
    'el_order': 2,
    'number_of_processors': 3,
    'current_control': 0,
    'Phi_vector': [1.5, None, 0.0, None],
    'freq': 130.0,
    'T': 60.0,
    't_step': 1.0,
    'phi': 0.0,
    'Signal_type': 'Rectangle',
    'Ampl_scale': 1.0,
    'CPE_activ': 1,
    'beta': 0.91,
    'K_A': 2621550.0,
    'beta_ground': 0.91,
    'K_A_ground': 0.0,
    'Full_Field_IFFT': 0,
    't_step_end': 1200,
    'VTA_from_divE': 0,
    'VTA_from_NEURON': 0,
    'VTA_from_E': 0,
    'Activation_threshold_VTA': 0,
    'spectrum_trunc_method': 'Octave Band Method',
    'trunc_param': 25000.0,
    'Truncate_the_obtained_full_solution': 0,
    'Show_paraview_screenshots': 0,

    'Solver_Type': 'BiCGSTAB',
    'FEniCS_MPI': 0,
    'Axon_Model_Type': 'McIntyre2002',
    'Approximating_Dimensions': [30.0, 30.0, 30.0],
    }

Results
--------

The results of both version ...

References
-----------

.. [Fonov2011] VS Fonov, AC Evans, K Botteron, CR Almli, RC McKinstry, DL Collins and BDCG, Unbiased average age-appropriate atlases for pediatric studies, NeuroImage,Volume 54, Issue 1, January 2011, ISSN 1053–8119, DOI: 10.1016/j.neuroimage.2010.07.033
