"""
Example of an ellipsoid.

We choose three radii:
    rx = 0.5
    ry = 0.25
    rz = 0.125

Now, we choose three different voxel numbers
in each direction: nx, ny, nz
Now, the voxels are anisotropic.

Outcome:
    The numpy array passed to NGSolve with
    shape (nx, ny, nz) works fine.

"""

import netgen.occ as occ
import numpy as np
from ngsolve import H1, Draw, GridFunction, Mesh, VoxelCoefficient, sqrt

Origin = (-1, -0.5, -0.25)
EndOfBox = (1, 0.5, 0.25)
edge_length = []
for i in range(len(Origin)):
    edge_length.append(EndOfBox[i] - Origin[i])

occgeo = occ.OCCGeometry(occ.Box(Origin, EndOfBox))
ngmesh = occgeo.GenerateMesh(maxh=0.1)
mesh = Mesh(ngmesh)
Draw(mesh)


rx = 0.5
ry = 0.25
rz = 0.125

nx = 88  # independent voxeldata mesh-size
ny = 88  # independent voxeldata mesh-size
nz = 88  # independent voxeldata mesh-size
data = np.zeros((nz, ny, nx))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            # x-, y-, z-position
            px, py, pz = (
                Origin[0] + i / (nx - 1) * edge_length[0],
                Origin[1] + j / (ny - 1) * edge_length[1],
                Origin[2] + k / (nz - 1) * edge_length[2],
            )
            # level set function of ellipsoid
            data[k][j][i] = 1.0 - sqrt(px**2 / rx**2 + py**2 / ry**2 + pz**2 / rz**2)

cf = VoxelCoefficient(Origin, EndOfBox, data, linear=False)

fes = H1(mesh, order=3)
gfu = GridFunction(fes)

gfu.Set(cf)

Draw(gfu)
