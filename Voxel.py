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

class DefaultMagneticResonanceImage():

    def __init__(self) -> None:
        self.__x = 50
        self.__y = 50
        self.__z = 50
        box = np.full((self.__x, self.__y, self.__z), 1)
        self.__array = np.empty((2*self.__x, 2*self.__y, 2*self.__z))
        self.__array[:self.__x, :self.__y, :self.__z] = box * 0
        self.__array[self.__x:, :self.__y, :self.__z] = box * 1
        self.__array[:self.__x, self.__y:, :self.__z] = box * 2
        self.__array[self.__x:, self.__y:, :self.__z] = box * 3
        self.__array[:self.__x, :self.__y, self.__z:] = box * 3
        self.__array[self.__x:, :self.__y, self.__z:] = box * 2
        self.__array[:self.__x, self.__y:, self.__z:] = box * 1
        self.__array[self.__x:, self.__y:, self.__z:] = box * 0

    def data_map(self) -> np.array:
        return self.__array

    def bounding_box(self) -> np.ndarray:
        starts = np.array([0, 0, 0])
        ends = starts + np.multiply(self._xyz_shape(), self._xyz_dimension())
        return tuple(np.array([starts, ends]))

    def _xyz_dimension(self) -> tuple:
        return (0.1, 0.1, 0.1)

    def _xyz_shape(self):
        return (2 * self.__x, 2 * self.__y, 2 * self.__z)


class MicroProbesCustomRodent():
    """MicroProbes Custom Rodent electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.

    Methods
    -------
    generate_geometry()
        Generate mesh of electrode.
    """

    # dimensions [mm]
    CONTACT_LENGTH = 0.01125
    LEAD_DIAMETER = 0.225
    TOTAL_LENGHTH = 13.3
    TUBE_THICKNESS = .01

    def __init__(self,
                 rotation: float = 0.0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        self.__translation = translation
        norm = np.linalg.norm(direction)
        self.__direction = tuple(direction / norm) if norm else (0, 0, 1)
        self.__boundaries = {'Body': 'Body', 'Contact_1': 'Contact_1'}

    def generate_geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        electrode = netgen.occ.Glue([self.__contact(), self.__body()])
        moved_electrode = electrode.Move(self.__translation)
        return moved_electrode

    def __body(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        tube_radius = radius + self.TUBE_THICKNESS
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=tube_radius)
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=tube_radius,
                                   h=self.TOTAL_LENGHTH - radius)
        point = tuple(np.array(self.__direction) * self.CONTACT_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        body = tip + lead - space
        body.bc(self.__boundaries['Body'])
        return body

    def __contact(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        radius = self.LEAD_DIAMETER * 0.5
        point = tuple(np.array(self.__direction) * self.CONTACT_LENGTH)
        space = netgen.occ.HalfSpace(p=point, n=self.__direction)
        center = tuple(np.array(self.__direction) * radius)
        tip = netgen.occ.Sphere(c=center, r=radius) * space
        lead = netgen.occ.Cylinder(p=center,
                                   d=self.__direction,
                                   r=radius,
                                   h=self.CONTACT_LENGTH - radius)

        if self.CONTACT_LENGTH <= radius:
            contact = tip
        else:
            contact = tip + lead
        contact.bc(self.__boundaries['Contact_1'])
        return contact

    def rename_boundaries(self, boundaries: dict) -> None:
        self.__boundaries.update(boundaries)

    def boundaries(self) -> dict:
        return self.__boundaries


class Ellipsoid:

    def __init__(self, start: tuple, end: tuple) -> None:
        self.__start = start
        self.__end = end

    def create(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        x, y, z = (np.array(self.__end) - np.array(self.__start)) / 2
        trasformator = netgen.occ.gp_GTrsf(mat=[x, 0, 0, 0, y, 0, 0, 0, z])
        sphere = netgen.occ.Sphere(c=netgen.occ.Pnt(1, 1, 1), r=1)
        ellipsoid = trasformator(sphere).Move(tuple(self.__start))
        ellipsoid.bc('Brain')
        return ellipsoid


mri = DefaultMagneticResonanceImage()
start = tuple(mri.bounding_box()[0])
end = tuple(mri.bounding_box()[1])
ellipsoid = Ellipsoid(start=start, end=end)
electrode = MicroProbesCustomRodent(translation=(5, 5, 5))
geometry = ellipsoid.create() - electrode.generate_geometry()
boundary_values = {'Body': 0, 'Contact_1': 1, }
mri_data = mri.data_map()
csf_pos = mri_data == 1
wm_pos = mri_data == 2
gm_pos = mri_data == 3
unknown_pos = mri_data == 0
data = np.zeros(mri.data_map().shape)
data[csf_pos] = 2.0
data[wm_pos] = 0.02
data[gm_pos] = 0.02
data[unknown_pos] = 0.02

mesh = netgen.occ.OCCGeometry(geometry).GenerateMesh()
ng_mesh = ngsolve.Mesh(mesh)
ng_mesh.Curve(order=2)
space = ngsolve.H1(mesh=ng_mesh,
                   order=2,
                   dirichlet='Brain|Contact_1|Body',
                   complex=False,
                   wb_withedges=False)

potential = ngsolve.GridFunction(space=space)
coefficient = ng_mesh.BoundaryCF(values=boundary_values)
potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)


# data = np.full((10, 10, 10), 0.2)
start = tuple(mri.bounding_box()[0])
end = tuple(mri.bounding_box()[1])
sigma = ngsolve.VoxelCoefficient(start=start,
                                 end=end,
                                 values=data,
                                 linear=False)

data = np.array([np.eye(3)]*25**3).reshape((25,25,25,3,3))
data = np.array([np.eye(3)])

diff = ngsolve.VoxelCoefficient(start=start,
                                 end=end,
                                 values=data,
                                 linear=False)

factor = diff * sigma

u = space.TrialFunction()
v = space.TestFunction()
a = ngsolve.BilinearForm(space=space, symmetric=True)
a += (factor * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx)
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

P = ngsolve.Integrate(ngsolve.grad(potential) *
                      ngsolve.Conj(sigma *
                                   ngsolve.grad(potential)),
                      ng_mesh)
print(1/P)