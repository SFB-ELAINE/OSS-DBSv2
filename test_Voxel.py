import numpy as np
import netgen
import ngsolve


# def voxel(x, y, z):
#     array = np.empty((x*2, y*2, z*2))

#     box = np.full((x, y, z), 1)
#     box1 = box
#     box2 = box * 2
#     box3 = box * 3
#     box4 = box * 4
#     box5 = box * 5
#     box6 = box * 6
#     box7 = box * 7
#     box8 = box * 8

#     array[:x, :y, :z] = box1
#     array[x:, :y, :z] = box2
#     array[:x, y:, :z] = box3
#     array[x:, y:, :z] = box4
#     array[:x, :y, z:] = box5
#     array[x:, :y, z:] = box6
#     array[:x, y:, z:] = box7
#     array[x:, y:, z:] = box8

#     return array

# x = 20
# y = 30
# z = 40
# start = (0, 0, 0)
# end = (2*x, 2*y, 2*z)
# model = netgen.occ.Box(p1=start, p2=end)
# geometry = netgen.occ.OCCGeometry(model)
# mesh = ngsolve.Mesh(geometry.GenerateMesh())
# data = voxel(x, y, z)
# cf = ngsolve.VoxelCoefficient(start=start, end=end, values=data, linear=True)
# ngsolve.Draw(cf, mesh, 'foo')

# # # Interpolation to zeroth order L2
# fes = ngsolve.L2(mesh, order=0, autoupdate=True)
# gfu = ngsolve.GridFunction(fes, autoupdate=True)
# gfu.Set(cf)

# ngsolve.Draw(gfu)

# refinements = 4
# for _ in range(refinements):
#     # mark all elements flagged with 1 for refinement
#     flags = np.isclose(gfu.vec.FV().NumPy(), 1.0)
#     for index, element in enumerate(mesh.Elements()):
#         mesh.SetRefinementFlag(ei=element, refine=flags[index])
#     mesh.Refine()
#     gfu.Set(cf)
#     ngsolve.Redraw()


geometry = netgen.occ.Sphere((0, 0, 0), 10)
geometry.bc('Box')
mesh = netgen.occ.OCCGeometry(geometry).GenerateMesh()
ng_mesh = ngsolve.Mesh(mesh)

space = ngsolve.H1(mesh=ng_mesh,
                   order=2,
                   dirichlet='Box',
                   complex=complex,
                   wb_withedges=False)

potential = ngsolve.GridFunction(space=space)
coefficient = ng_mesh.BoundaryCF(values={'Box': 1})
potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)

data = np.full((10, 10, 10), 0.2)
sigma = ngsolve.VoxelCoefficient(start=(0, 0, 0),
                                 end=(10, 10, 10),
                                 values=data,
                                 linear=False)

data = np.array([np.eye(3)] * 25**3).reshape((25,25,25,3,3))
#data = np.array([np.eye(3)])
diff = ngsolve.VoxelCoefficient(start=(0, 0, 0),
                                 end=(10, 10, 10),
                                 values=data,
                                 linear=False)
factor = diff * sigma

u = space.TrialFunction()
v = space.TestFunction()
a = ngsolve.BilinearForm(space=space, symmetric=True)
a += (diff * sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx)
f = ngsolve.LinearForm(space=space)
preconditioner = ngsolve.Preconditioner(bf=a, type="bddc", coarsetype="h1amg")


a.Assemble()
f.Assemble()
inverse = ngsolve.CGSolver(mat=a.mat,
                           pre=preconditioner.mat,
                           printrates=True,
                           maxsteps=10000,
                           precision=1e-12)
r = f.vec.CreateVector()
r.data = f.vec - a.mat * potential.vec
potential.vec.data = potential.vec.data + inverse * r
