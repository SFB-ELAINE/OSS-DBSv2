import json
from pprint import pprint

import ngsolve

import ossdbs

ngsolve.SetNumThreads(8)
input_file = "input_surface_impedance_lempka2009.json"
with open(input_file) as fp:
    input_dict = json.load(fp)
model_geometry = ossdbs.generate_model_geometry(input_dict)
with ngsolve.TaskManager():
    mesh = ossdbs.generate_mesh(input_dict)
    mesh.curve(2)
    boundary_areas = mesh.get_boundary_areas()
pprint(boundary_areas)
print("#####")
pprint(model_geometry.contacts)
