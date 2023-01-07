import ngsolve
import netgen
from netgen.occ import Sphere, Pnt, Vec
import json
import jsons
from dataclasses import dataclass
from tests.models_data import Vertex_data,CurveRepresentation,Edge_data,Wire_data,Face_data,Solid_data
from ossdbs.electrodes.abbott_stjude_active_tip_6142_6145 import AbbottStjudeActiveTip6142_6145
from tests.static_classes import EnhancedJSONEncoder

import dataclasses, json


sp1 = Sphere(Pnt(0,0,0), 0.2)
sp2 = Sphere(Pnt(0.5,0,0), 0.2)
sp3 = Sphere(Pnt(0,0,0.5), 0.2)
sp4 = Sphere(Pnt(0,0.2,0.7), 0.2)

vector = Vec(1, 1, 1)
all = sp1 + sp2 + sp3 + sp4
compound = all.Move(vector)


def convert_vertex(vertex) -> tuple:
    vertex_data = json.loads("{"+str(vertex)+"}")
    transformation = vertex_data['Location']['Transformation']
    coordinates = vertex_data['TShape']['Pnt']['gp_Pnt']
    return Vertex_data(transformation,coordinates)

def convert_edge(edge) -> Edge_data:
    edge_data = json.loads("{"+str(edge)+"}")
    tshape = edge_data['TShape']
    curve = tshape['CurveRepresentation']
    new_curve=CurveRepresentation(curve['className'],
                                  curve['First'],
                                  curve['Last'],
                                  curve['UV1'],
                                  curve['UV2'],
                                  curve['PCurve'],
                                  curve['Surface'],
                                  curve['Location']['Transformation']['Matrix']
                                  )
    result=Edge_data(tshape['Flags'],
                    tshape['Orientable'],
                    tshape['Tolerance'],
                    new_curve)
    return result

def convert_wire(w) -> Wire_data:
    wire_data=json.loads('{'+str(w)+'}')
    wire_obj = Wire_data(wire_data['className'],
                wire_data['TShape'],
                wire_data['Location'],
                wire_data['Orient'],
                )
    return wire_obj

def convert_face(f)-> Face_data:
    face_obj=json.loads('{'+str(f)+'}')
    face_data=Face_data(face_obj['Orient'],
                        face_obj['Location']['Transformation']['Matrix'],
                        face_obj["TShape"]['Surface'])
    return face_obj

def convert_solid(solid_obj) -> Solid_data:
    solid=json.loads('{'+str(solid_obj)+'}')
    nbChildren = solid['TShape']['NbChildren']
    children = [convert_solid(i) for i in solid_obj.solids] if nbChildren>1 else []
    faces = [convert_face(f) for f in solid_obj.faces] if nbChildren==1 else []
    vertices = [convert_vertex(v) for v in solid_obj.vertices] if nbChildren==1 else []
    edges = [convert_edge(e) for e in solid_obj.edges] if nbChildren==1 else []
    wires = [convert_wire(w) for w in solid_obj.wires] if nbChildren==1 else []
    solid_data = Solid_data(children,
                            solid['Orient'],
                            nbChildren, 
                            solid['Location']['Transformation']['Matrix'],
                            faces,
                            vertices,
                            edges,
                            wires
                            )
    return solid_data

def solid_to_json(s):
    return json.dumps(convert_solid(s),cls=EnhancedJSONEncoder)

assert([solid_to_json(all)==solid_to_json(i) for i in all.solids]==[False]*len(all.solids))
assert(solid_to_json(all)!=solid_to_json(compound))

def generate_json_AbbottStjudeActiveTip6142_6145():
    TESTDATA = [
        ((0,  0, 0), (0, 0, 1), 0.0), # UseCase 1,2,3
        ((1, -2, 3), (0, 0, 1), 0.0), # UseCase 4,5,6
        ((1, -2, 3), (2, 0, 1), 0.0), # UseCase 7,8,9
    ]
    for i in range(len(test_cases)):
        translation,direction,rotation=test_cases[i]
        filename="AbbottStjudeActiveTip6142_6145_"+str(i)+".json"
        electrode = AbbottStjudeActiveTip6142_6145(translation,
                                                   direction,
                                                   rotation)
        with open(filename, 'w') as outfile:
            json.dump(convert_solid(electrode), outfile, cls=EnhancedJSONEncoder)

def main():
    generate_json_AbbottStjudeActiveTip6142_6145()


if __name__ == '__main__':
    main()
