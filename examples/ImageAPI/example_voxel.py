"""
Example of an ellipsoid.

We choose three radii:
    rx = 0.5
    ry = 0.25
    rz = 0.125

Outcome:
    The numpy array passed to NGSolve with
    shape (n, n, n), n being the number of voxels,
    is filled as data[z][y][x]

"""

import netgen.occ as occ
import numpy as np
from ngsolve import H1, Draw, GridFunction, Mesh, VoxelCoefficient, sqrt

occgeo = occ.OCCGeometry(occ.Box((-1, -1, -1), (1, 1, 1)))
ngmesh = occgeo.GenerateMesh(maxh=0.3)
mesh = Mesh(ngmesh)
Draw(mesh)

Origin = (-1, -1, -1)
edge_length = (2, 2, 2)

rx = 0.5
ry = 0.25
rz = 0.125

n = 22  # independent voxeldata mesh-size
data = np.zeros((n, n, n))
for i in range(n):
    for j in range(n):
        for k in range(n):
            # x-, y-, z-position
            px, py, pz = (
                Origin[0] + i / (n - 1) * edge_length[0],
                Origin[1] + j / (n - 1) * edge_length[1],
                Origin[2] + k / (n - 1) * edge_length[2],
            )
            # level set function of ellipsoid
            data[i][j][k] = 1.0 - sqrt(px**2 / rx**2 + py**2 / ry**2 + pz**2 / rz**2)

cf = VoxelCoefficient((-1, -1, -1), (1, 1, 1), data, linear=False)

fes = H1(mesh, order=2)
gfu = GridFunction(fes)

gfu.Set(cf)

Draw(gfu)
