"""
Example of an MRI image.

The MRI image is loaded, a VoxelCoefficient
function is created and voxel space.
The affine transformation is used to map
real space to voxel space.

Outcome:
    The MRI image passed to NGSolve
    is reproduced correctly.

"""

import netgen.occ as occ
import nibabel
import numpy as np
from ngsolve import (
    H1,
    CoefficientFunction,
    Draw,
    GridFunction,
    Mesh,
    VoxelCoefficient,
    x,
    y,
    z,
)
from numpy.linalg import inv

#############################
# PROCESS IMAGE

# Example MRI image loaded
mri_image = nibabel.load("segmask.nii.gz")
# Get the dimensions of the voxel data
nx, ny, nz = mri_image.header.get_data_shape()
affine = mri_image.affine
trafo_matrix = np.array([affine[0, :3], affine[1, :3], affine[2, :3]])
translation = affine[:3, 3]

end_of_box = trafo_matrix.dot([nx, ny, nz]) + translation
print(end_of_box)
# load data
data = mri_image.get_fdata()
# swap x and z axes for VoxelCoefficient
data = np.swapaxes(data, 0, 2)

Origin = tuple(translation)
EndOfBox = tuple(end_of_box)

OriginVoxelSpace = (0, 0, 0)
EndOfBoxVoxelSpace = (nx, ny, nz)

occgeo = occ.OCCGeometry(occ.Box(Origin, EndOfBox))
ngmesh = occgeo.GenerateMesh(maxh=10.0)
mesh = Mesh(ngmesh)
Draw(mesh)


fes = H1(mesh, order=3)
gfu = GridFunction(fes)
# Extract trafo_matrix from header
# ravel because ngsolve needs a vector in a tuple
inv_trafo_matrix = tuple(inv(trafo_matrix).ravel())
mm_space_coordinates = CoefficientFunction((x, y, z), dims=(3,))
mm_space_to_voxel_space = CoefficientFunction(inv_trafo_matrix, dims=(3, 3))
# Getting image offset for translation
translation = CoefficientFunction(tuple(translation), dims=(3,))
trafocf = mm_space_to_voxel_space * (mm_space_coordinates - translation)
cf = VoxelCoefficient(
    OriginVoxelSpace, EndOfBoxVoxelSpace, data, linear=False, trafocf=trafocf
)
gfu.Set(cf)
Draw(gfu, mesh, "MRI")
