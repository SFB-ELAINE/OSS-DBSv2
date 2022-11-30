import numpy as np
import netgen
import ngsolve


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


# Create Geometry and Mesh

mri = DefaultMagneticResonanceImage()
start = tuple(mri.bounding_box()[0])
end = tuple(mri.bounding_box()[1])

ellipsoid = Ellipsoid(start=start, end=end)
electrode = MicroProbesCustomRodent(translation=(5, 5, 5))
electrode.rename_boundaries({'Body': 'E1B', 'Contact_1': 'E1C1'})
electrode2 = MicroProbesCustomRodent(translation=(4.5, 4.5, 5))
electrode.rename_boundaries({'Body': 'E2B', 'Contact_1': 'E2C1'})
geometry = ellipsoid.create() - electrode.generate_geometry()
geometry = geometry - electrode2.generate_geometry()

mesh = netgen.occ.OCCGeometry(geometry).GenerateMesh()
ng_mesh = ngsolve.Mesh(mesh)
ng_mesh.Curve(order=2)


# Create VoxelCoefficient-Function for Conductivity

mri_data = mri.data_map()
csf_pos = mri_data == 1
wm_pos = mri_data == 2
gm_pos = mri_data == 3

mri_data = np.full(mri.data_map().shape, 0.02)
mri_data[csf_pos] = 2.0
mri_data[wm_pos] = 0.02
mri_data[gm_pos] = 0.02

shape = (*mri_data.shape, 3, 3)

start = tuple(mri.bounding_box()[0])
end = tuple(mri.bounding_box()[1])
sigma = ngsolve.VoxelCoefficient(start=start,
                                 end=end,
                                 values=mri_data,
                                 linear=False)

diff = ngsolve.VoxelCoefficient(start=start,
                                end=end,
                                values=np.array([np.eye(3)]),
                                linear=False)
sigma = sigma * diff

mri_data = np.array([np.eye(3) * d for d in mri_data.flatten()]).reshape(shape)
# sigma = ngsolve.VoxelCoefficient(start=start,
#                                  end=end,
#                                  values=mri_data,
#                                  linear=False)


# Refine Mesh

grid_function = ngsolve.GridFunction(space=ngsolve.L2(ng_mesh, order=0))
cf = ngsolve.VoxelCoefficient(start=tuple(start),
                              end=tuple(end),
                              values=csf_pos.astype(float),
                              linear=False)
grid_function.Set(cf)
flags = grid_function.vec.FV().NumPy()
for index, element in enumerate(ng_mesh.Elements()):
    ng_mesh.SetRefinementFlag(ei=element, refine=flags[index])

ng_mesh.Refine()
ng_mesh.Curve(2)

space = ngsolve.H1(mesh=ng_mesh,
                   order=2,
                   dirichlet='E1B|E2B|E1C1|E2C1',
                   complex=False,
                   wb_withedges=False)


# Solve Problem

boundary_values = {'E1B': 0, 'E1C1': 1, 'E2B': 0, 'E2C1': 1, }
potential = ngsolve.GridFunction(space=space)
coefficient = ng_mesh.BoundaryCF(values=boundary_values)
potential.Set(coefficient=coefficient, VOL_or_BND=ngsolve.BND)

#with ngsolve.TaskManager():

u = space.TrialFunction()
v = space.TestFunction()
a = ngsolve.BilinearForm(space=space, symmetric=True)
a += (sigma * ngsolve.grad(u) * ngsolve.grad(v) * ngsolve.dx)
f = ngsolve.LinearForm(space=space)
preconditioner = ngsolve.Preconditioner(bf=a, type="bddc", coarsetype="local")
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
                                   ngsolve.grad(potential)), ng_mesh)
print(1/P)