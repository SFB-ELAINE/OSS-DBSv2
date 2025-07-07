"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and
local mesh refinements are taken
into account during the meshing
process.
"""

import json

import ngsolve

import ossdbs
from ossdbs.utils.vtk_export import FieldSolution

ossdbs.set_logger()

# enter file name as needed
filename = "oss-dbs_parameters.json"


def save(cf, filename: str) -> None:
    """Save solution in VTK format."""
    names = [f"{cf.label}_real"]
    if cf.is_complex:
        names.append(f"{cf.label}_imag")

    coefficients = [cf.solution.real]
    if cf.is_complex:
        coefficients.append(cf.solution.imag)

    vtk = ngsolve.VTKOutput(
        ma=cf.mesh,
        coefs=coefficients,
        names=names,
        filename=filename,
        subdivision=0,
    )
    vtk.Do(vb=ngsolve.BND)


with open(filename) as fp:
    settings = json.load(fp)

electrodes = ossdbs.generate_electrodes(settings)
electrode = electrodes[0]
electrode.export_electrode(
    output_path=".", brain_dict=settings["BrainRegion"], n_electrode=0
)

# generate mesh with brain
mesh = ossdbs.generate_mesh(settings)
mesh.save("test.vol.gz")

print(mesh.boundaries)
print(mesh.materials)

bnd_dict = {}
for idx in range(1, 9):
    bnd_dict[f"E1C{idx}"] = idx

cf = FieldSolution(
    solution=mesh.boundary_coefficients(bnd_dict),
    label="bnd",
    mesh=mesh.ngsolvemesh,
    is_complex=False,
)
save(cf, "electrode_in_brain")
