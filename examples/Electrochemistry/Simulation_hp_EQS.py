import impedancefitter as ifit
import netgen.occ as occ
import numpy as np
import pandas as pd
from ngsolve import (
    BND,
    H1,
    VOL,
    BilinearForm,
    BoundaryFromVolumeCF,
    CGSolver,
    CoefficientFunction,
    Conj,
    Draw,
    GridFunction,
    HDiv,
    Integrate,
    LinearForm,
    Mesh,
    Preconditioner,
    TaskManager,
    ds,
    dx,
    grad,
    specialcf,
    sqrt,
)
from ngsolve.utils import Laplace

import ossdbs

ecm = ifit.get_equivalent_circuit_model("CPE")
Z = ecm.eval(omega=2.0 * np.pi * 1e3, alpha=0.5, k=1e5)  # impedance in Ohm
sigma = 0.1 * 1e-3  # conductivity in S/mm
Vlow = 0
Vhigh = 1.0
contact_high = "Contact_1"
contact_low = "Contact_8"

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

cfvals = {contact_low: Vlow, contact_high: Vhigh}
cf = mesh.BoundaryCF(cfvals, default=-1)
Draw(cf, mesh, "BND")
area1 = Integrate(CoefficientFunction(1.0) * ds(contact_high), mesh)
ys = 1.0 / (Z * area1)
print("Surface admittance [1 / (Ohm * mm^2)]: ", ys)
print("Surface resistance [Ohm * mm^2]: ", 1.0 / ys)

fes = H1(
    mesh,
    order=order,
    dirichlet=contact_low,
    autoupdate=True,
    wb_withedges=False,
    complex=True,
)

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes, symmetric=True)
a += Laplace(coef=sigma)
# Robin BC
a += CoefficientFunction(ys) * u * v * ds(contact_high)

f = LinearForm(fes)
f += (
    CoefficientFunction(ys)
    * CoefficientFunction(cfvals[contact_high])
    * v
    * ds(contact_high)
)

c = Preconditioner(a, type="bddc", coarsetype="h1amg")

gfu = GridFunction(fes, autoupdate=True)

# dict for results
results_dict = {}
results_dict["I"] = []
results_dict["I1"] = []
results_dict["I2"] = []
results_dict["dof"] = []
results_dict["nels"] = []
results_dict["Volume"] = []
results_dict["Field"] = []
results_dict["maxerr"] = []
results_dict["errorestim"] = []

space_flux = HDiv(mesh, order=order - 1, autoupdate=True, complex=True)
gf_flux = GridFunction(space_flux, "flux", autoupdate=True)


def CalcError():
    """Estimate error by ZZ error estimator."""
    flux = grad(gfu)
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
    gfu.Set(cf, BND)
    res = gfu.vec.CreateVector()
    res.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += inv * res

    # use that U = 1, P = U^2 / R, I = U / R = P / U = P (if U = 1 as here)
    # Integrate to compute P
    I_vol = 1e3 * sigma * Integrate(-grad(gfu) * Conj(-grad(gfu)), mesh)
    print("I = ", I_vol)
    tmp = 1e3 * Integrate(
        ys * (cf - gfu) * Conj(cf - gfu),
        mesh,
        definedon=mesh.Boundaries(contact_high),
    )
    print("tmp = ", tmp)
    I_vol += tmp
    # account for voltage drop
    I_vol /= Vhigh - Vlow

    n = specialcf.normal(3)
    Jn = BoundaryFromVolumeCF(-grad(gfu)) * BoundaryFromVolumeCF(n)
    I1 = 1e3 * sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_high))
    I2 = 1e3 * sigma * Integrate(Jn, mesh, definedon=mesh.Boundaries(contact_low))
    # evalaute field at bottom of well
    mip = mesh(0.0, 0.0, 1e-7)
    field = sqrt(grad(gfu) * grad(gfu))
    field = field(mip) * 1e3
    volume = Integrate(CoefficientFunction(1.0) * dx, mesh)
    errorestim, maxerr = CalcError()
    print("current [mA] = ", np.round(I_vol, 4))
    print("Impedance [Ohm] = ", np.round(1e3 / I_vol, 4))
    print("current 1 [mA] = ", I1)
    print("current 2 [mA] = ", I2)
    print("field [V/m] = ", field)
    results_dict["dof"].append(fes.ndof)
    results_dict["I"].append(I_vol)
    results_dict["I1"].append(I1)
    results_dict["I2"].append(I2)
    results_dict["Volume"].append(volume)
    results_dict["Field"].append(field)
    results_dict["nels"].append(mesh.ne)
    results_dict["errorestim"].append(errorestim)
    results_dict["maxerr"].append(maxerr)

Draw(gfu)
field = -grad(gfu)
Draw(field, mesh, "field")

data = pd.DataFrame(results_dict)
data.to_csv("results_Z_EQS.csv", index=False)

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
