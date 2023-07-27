"""
A single micro probe rodent electrode 
with an encapsulation layer is generated using custom parameters
"""

import ossdbs
from ngsolve import Draw, Mesh, TaskManager
import netgen.occ as occ

parameters ={'exposed': 0.15,
             'contact_radius': 0.1125,
              'lead_radius':  0.1125,
              'total_length': 13.3, 
                'wire_radius' :0.09 }


settings = {"Electrodes":
            [{"Name": "MicroProbesRodentElectrodeCustom",
              'CustomParameters': parameters,
              "Rotation[Degrees]": 0,
              "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
              "TipPosition": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 0},
              "EncapsulationLayer":
                {"Thickness[mm]": 0.1,
                 }
              },
             ]
            }

electrodes = ossdbs.generate_electrodes(settings)
electrode_settings = settings["Electrodes"][0]
electrode = electrodes[0]
encap = electrode.encapsulation_geometry(electrode_settings["EncapsulationLayer"]["Thickness[mm]"])
occgeo = occ.OCCGeometry(occ.Glue([electrode.geometry, encap]))
with TaskManager():
    mesh = Mesh(occgeo.GenerateMesh())
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())

# Clip in the center and show the scalar function to see the layer
Draw(mesh.MaterialCF({"EncapsulationLayer": 1.0}), mesh, "Material")
