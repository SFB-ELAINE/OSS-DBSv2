Image API examples
==================

[example_voxel.py](examples/ImageAPI/example_voxel.py)
------------------------------------------------------

In this initial example, we can visualize the issue with how data in the Nifti image is stored and how NGSolve reads in the data. We see that `data` is stored in the Nifti file format in the order of `(x1,y1,z1), (x1,y1,z2),...(x1,y2,z1)`. The issue lies that the data is being read by NGSolve in the `(x1,y1,z1),(x2,y1,z1),...(x1,y2,z1)` order. Here, we use an ellipsoid as an example to visualize how the x-axis and z-axis are flipped. 



[example_voxel_nx_ny_nz.py](examples/ImageAPI/example_voxel_nx_ny_nz.py)
------------------------------------------------------------------------

In this example, we fix the aformentioned issue by reshaping and storing the data in the order that NGSolve would traverse the data in: `(x1,y1,z1),(x2,y1,z1),...(x1,y2,z1)`. Here, we also use an ellipsoid as an example to visualize how the x-axis and z-axis are not flipped. Note, that this example also works when nz, nx, and ny are different dimensions.

[example_voxel_nx_ny_nz_anisotropic.py](examples/ImageAPI/example_voxel_nx_ny_nz_anisotropic.py)
------------------------------------------------------------------------------------------------

In this example, we use anisotropic voxels to test the VoxelCoefficient function. The results show that this function works even for non-cubicle voxel sizes.

[example_voxel_nx_ny_nz_trafocf.py](examples/ImageAPI/example_voxel_nx_ny_nz_trafocf.py)
----------------------------------------------------------------------------------------

In this example, we test the `trafocf` input parameter by creating a `trafo_matrix` that is used to transform the points in mm space to voxel space. The `trafocf` input is then passed into the VoxelCoefficient function so points in mm space can be mapped to `data` in voxel space.

We can clearly visualize this transformation in the z-y plane after using the trafocf parameter with the inverted matrix. Even though the orginal ry and rz values differ, we can see that after applying the transformation, the y and z radius of the ellipsoid are the same.


Test DTI images
===============

Download for example the Lead-DBS data (link on their website or the README of their github repository).
We used the following DTI image: `templates/space/MNI_ICBM_2009b_NLIN_ASYM/IITmean_tensor_Norm_mapping.nii.gz`


