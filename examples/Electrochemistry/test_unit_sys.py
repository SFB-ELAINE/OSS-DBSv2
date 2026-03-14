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

Z = 1e1
sigma = 0.1  # conductivity in S/mm
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
finalgeo = finalgeo.Scale((0, 0, 0), 1e-3)
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

fes_list = []
for _ in range(2):
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
for i in range(2):
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
a += sum_u * v[-1] * dx
a += sum_v * u[-1] * dx

f = LinearForm(fes)

current = 0.002
for i, electrode in enumerate(electrodes):
    length = Integrate(CoefficientFunction(1) * ds(electrode), mesh)
    print(current / length)
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
    gfu.components[0].Set(0)
    inv = CGSolver(mat=a.mat, pre=c.mat, printrates=True, maxsteps=100, precision=1e-8)
    res = f.vec.CreateVector()
    res.data = f.vec
    gfu.vec.data = gfu.vec.data + inv * f.vec

    n = specialcf.normal(3)
    Jn = BoundaryFromVolumeCF(-grad(gfu.components[0])) * BoundaryFromVolumeCF(n)
    I1 = sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_high))
    I2 = sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_low))
    # evalaute field at bottom of well
    mip = mesh(0.0, 15.0e-3, 25e-3)
    field = sqrt(grad(gfu.components[0]) * grad(gfu.components[0]))
    field = field(mip) * 1e3
    volume = Integrate(CoefficientFunction(1.0) * dx, mesh)
    errorestim, maxerr = CalcError()
    print("current 1 [mA] = ", I1)
    print("current 2 [mA] = ", I2)
    print("field [V/m] = ", field)
    results_dict["dof"].append(fes.ndof)
    results_dict["I1"].append(I1)
    results_dict["I2"].append(I2)
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
for i in range(2):
    print(gfu.components[i + 1].vec)

"""
vtk = VTKOutput(
    ma=mesh,
    coefs=[gfu, -grad(gfu) * 1e3],
    names=["potential", "field"],
    filename="result_manual",
    subdivision=0,
)
# Exporting the results:
vtk.Do()
"""
