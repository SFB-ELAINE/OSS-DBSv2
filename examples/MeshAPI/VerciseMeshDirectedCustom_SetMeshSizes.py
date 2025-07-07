"""
Example of a full model geometry
with one Vercise electrode without
encapsulation layer.
The mesh is generated and
local mesh refinements are taken
into account during the meshing
process.
"""

import ngsolve

import ossdbs
from ossdbs.utils.vtk_export import FieldSolution


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


settings = {
    "Electrodes": [
        {
            "Name": "BostonScientificVerciseDirectedCustom",
            "CustomParameters": {
                "tip_length": 1.5,
                "contact_length": 0.5,
                "contact_spacing": 0.5,
                "lead_diameter": 1.3,
                "total_length": 450.0,
            },
            "Rotation[Degrees]": 0,
            # "Direction": {"x[mm]": -0.5775, "y[mm]": -0.5775, "z[mm]": -0.5775},
            # could not be built?
            "Direction": {"x[mm]": 0.7075, "y[mm]": 0.0, "z[mm]": 0.7075},
            # "Direction": {"x[mm]": 1, "y[mm]": 1, "z[mm]": 1},
            # "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
            "TipPosition": {"x[mm]": 0, "y[mm]": 15, "z[mm]": -3},
            "EncapsulationLayer": {
                "Thickness[mm]": 0.0,  # indicates that no encapsulation is modelled
            },
            "Contacts": [
                {
                    "Contact_ID": 1,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 2,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 3,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 4,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 5,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 6,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 7,
                    "MaxMeshSizeEdge": 0.01,
                },
                {
                    "Contact_ID": 8,
                    "MaxMeshSizeEdge": 0.01,
                },
            ],
        }
    ],
    "BrainRegion": {
        "Center": {"x[mm]": 5, "y[mm]": 14, "z[mm]": -4.5},
        "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
        "Shape": "Ellipsoid",
    },
    "Mesh": {"LoadMesh": False, "SaveMesh": False},
    "ExportElectrode": False,
}

# generate mesh with brain
with ngsolve.TaskManager():
    mesh = ossdbs.generate_mesh(settings)

print(mesh.boundaries)
print(mesh.materials)

bnd_dict = {}
for idx in range(1, 9):
    bnd_dict[f"E1C{idx}"] = idx

ngsolve.Draw(mesh.boundary_coefficients(bnd_dict), mesh.ngsolvemesh, "bnd")

cf = FieldSolution(
    solution=mesh.boundary_coefficients(bnd_dict),
    label="bnd",
    mesh=mesh.ngsolvemesh,
    is_complex=False,
)
save(cf, "electrode_custom_in_brain")
