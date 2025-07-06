"""
Example of reading in a DTI image and using get_component function for xy and zy
Tests function for randomly generated 100 voxel points.

Outcome:
    get_component returns 3D array for each component with correct values
"""

import random

from ossdbs.utils.nifti1image import DiffusionTensorImage

#############################


# Loading DTI image
dti_image = DiffusionTensorImage("IITMeanTensor_NormMapping.nii.gz")
dti_data_all = dti_image.data
diffusion_xx = dti_image.get_component("xx")
diffusion_xy = dti_image.get_component("xy")
diffusion_xz = dti_image.get_component("xz")
diffusion_yx = dti_image.get_component("yx")
diffusion_yy = dti_image.get_component("yy")
diffusion_yz = dti_image.get_component("yz")
diffusion_zx = dti_image.get_component("zx")
diffusion_zy = dti_image.get_component("zy")
diffusion_zz = dti_image.get_component("zz")

dti_shape = dti_image.xyz_shape
print("Dimensions of component data:", diffusion_xy.shape)
# Testing by referencing specific voxels
# Randomly generate integer indices in range x,y,z baed on shape of data
test_indices = [
    (
        random.randint(0, dti_shape[0] - 1),
        random.randint(0, dti_shape[1] - 1),
        random.randint(0, dti_shape[2] - 1),
    )
    for i in range(25)
]
for voxel in test_indices:
    print("Test Voxel:", voxel)
    print("Value from get_component (xy):", diffusion_xy[voxel])
    print("Value from indexing data map:", dti_data_all[voxel][1])
    print("Value from get_component (zy):", diffusion_zy[voxel])
    print("Value from indexing data map:", dti_data_all[voxel][4], "\n")
