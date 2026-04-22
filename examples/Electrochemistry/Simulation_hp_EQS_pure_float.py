import impedancefitter as ifit
import netgen.occ as occ
import numpy as np
import pandas as pd
from ngsolve import (
    H1,
    VOL,
    BilinearForm,
    BoundaryFromVolumeCF,
    CGSolver,
    CoefficientFunction,
    Conj,
    Draw,
    FESpace,
    GridFunction,
    HDiv,
    Integrate,
    LinearForm,
    Mesh,
    NumberSpace,
    Preconditioner,
    TaskManager,
    ds,
    dx,
    grad,
    specialcf,
    sqrt,
)

import ossdbs

ecm = ifit.get_equivalent_circuit_model("R")
Z = ecm.eval(omega=2.0 * np.pi * 1e3, R=1e3)  # impedance in Ohm
sigma = 0.1 * 1e-3  # conductivity in S/mm
contact_high = "Contact_1"
contact_low = "Contact_8"
electrodes = [contact_high, contact_low]
# FEM order
order = 2

brain_geo = occ.Sphere(c=(0, 0, 0), r=50)
vercise = ossdbs.electrodes.defaults.BostonScientificVercise()
finalgeo = brain_geo - vercise.geometry
for e in finalgeo.edges:
    if e.name is not None:
        if contact_high in e.name or contact_low in e.name:
            e.hpref = 1

occgeo = occ.OCCGeometry(finalgeo)
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
    mesh.RefineHP(levels=2)
    mesh.Curve(order)
Draw(mesh)

area1 = Integrate(CoefficientFunction(1.0) * ds(contact_high), mesh)
ys = 1.0 / (np.real(Z) * area1)
print("Surface admittance [1 / (Ohm * mm^2)]: ", ys)
print("Surface resistance [Ohm * mm^2]: ", 1.0 / ys)

h1 = H1(
    mesh,
    order=order,
    autoupdate=True,
    wb_withedges=False,
    complex=False,
)

n_float = len(electrodes)
fes_list = []
for _ in range(n_float):
    fes_list.append(
        NumberSpace(
            mesh,
            complex=False,
        )
    )
# Lagrange multiplier to enforce sum of floating potentials = 0
fes_list.append(
    NumberSpace(
        mesh,
        complex=False,
    )
)

fes = FESpace([h1, *fes_list])

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += sigma * grad(u[0]) * grad(v[0]) * dx

sum_u = None
sum_v = None
for i in range(n_float):
    # Robin BC
    a += (
        CoefficientFunction(ys)
        * (u[0] - u[i + 1])
        * (v[0] - v[i + 1])
        * ds(electrodes[i])
    )
    if i == 0:
        sum_u = u[1]
        sum_v = v[1]
    else:
        sum_u += u[i + 1]
        sum_v += v[i + 1]

# Lagrange multiplier constraint: sum of floating potentials = 0
lam = u[1 + n_float]
mu = v[1 + n_float]
bnd = electrodes[0]
a += sum_u * mu * ds(bnd)
a += sum_v * lam * ds(bnd)
# Small regularization for iterative solver
eps = 1e-12
a += eps * lam * mu * ds(bnd)

f = LinearForm(fes)

current = 0.002
for i, electrode in enumerate(electrodes):
    length = Integrate(CoefficientFunction(1.0) * ds(electrode), mesh)
    f += (-1) ** i * current / length * v[i + 1] * ds(electrode)

c = Preconditioner(a, type="local")

gfu = GridFunction(fes, autoupdate=True)

# dict for results
results_dict = {}
results_dict["I1"] = []
results_dict["I2"] = []
results_dict["dof"] = []
results_dict["nels"] = []
results_dict["Volume"] = []
results_dict["Field"] = []
results_dict["maxerr"] = []
results_dict["errorestim"] = []
results_dict["V1"] = []
results_dict["V2"] = []

space_flux = HDiv(mesh, order=order - 1, autoupdate=True, complex=False)
gf_flux = GridFunction(space_flux, "flux", autoupdate=True)


def CalcError():
    """Estimate error by ZZ error estimator."""
    flux = grad(gfu.components[0])
    # interpolate finite element flux into H(div) space:
    gf_flux.Set(flux)

    # Gradient-recovery error estimator
    err = (flux - gf_flux) * Conj(flux - gf_flux)
    elerr = Integrate(err.real, mesh, VOL, element_wise=True)

    maxerr = max(elerr)
    print("maxerr = ", maxerr)

    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, elerr[el.nr] > 0.5 * maxerr)
    return sqrt(sum(elerr)), maxerr


with TaskManager():
    print("Assemble a")
    a.Assemble()
    print("Assemble b")
    f.Assemble()
    inv = CGSolver(
        mat=a.mat, pre=c.mat, printrates=True, maxsteps=10000, precision=1e-12
    )
    res = gfu.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += inv * res

    n = specialcf.normal(3)
    Jn = BoundaryFromVolumeCF(-grad(gfu.components[0])) * BoundaryFromVolumeCF(n)
    I1 = 1e3 * sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_high))
    I2 = 1e3 * sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_low))
    # floating potentials
    V1 = gfu.components[1].vec[0]
    V2 = gfu.components[2].vec[0]
    # evaluate field at bottom of well
    mip = mesh(0.0, 0.0, 1e-7)
    field = sqrt(grad(gfu.components[0]) * grad(gfu.components[0]))
    field = field(mip) * 1e3
    volume = Integrate(CoefficientFunction(1.0) * dx, mesh)
    errorestim, maxerr = CalcError()
    print("current 1 [mA] = ", I1)
    print("current 2 [mA] = ", I2)
    print("floating potential 1 [V] = ", V1)
    print("floating potential 2 [V] = ", V2)
    print("field [V/m] = ", field)
    results_dict["dof"].append(fes.ndof)
    results_dict["I1"].append(I1)
    results_dict["I2"].append(I2)
    results_dict["V1"].append(V1)
    results_dict["V2"].append(V2)
    results_dict["Volume"].append(volume)
    results_dict["Field"].append(field)
    results_dict["nels"].append(mesh.ne)
    results_dict["errorestim"].append(errorestim)
    results_dict["maxerr"].append(maxerr)

Draw(gfu.components[0])
field = -grad(gfu.components[0])
Draw(field, mesh, "field")

data = pd.DataFrame(results_dict)
data.to_csv("results_Z_pure_float.csv", index=False)
