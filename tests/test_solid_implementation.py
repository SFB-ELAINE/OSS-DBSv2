from netgen.occ import *
import json

sp1 = Sphere (Pnt(0,0,0), 0.2)
sp2 = Sphere (Pnt(0.5,0,0), 0.2)
sp3 = Sphere (Pnt(0,0,0.5), 0.2)
sp4 = Sphere (Pnt(0,0.2,0.7), 0.2)
#
v=Vec(1,1,1)
all = sp1+sp2+sp3+sp4
b=all.Move(v)

def convert_vertex(v):
    return ((s:=json.loads("{"+str(v)+"}"))['Location']['Transformation'],s['TShape']['Pnt']['gp_Pnt']) 
def convert_edge(e):
    obj=json.loads("{"+str(e)+"}")
    orient=obj['Orient']
    location=obj['Location']
    tshape=obj['TShape']
    curve=tshape['CurveRepresentation']
    curve['Matrix']=curve['Location']['Transformation']['Matrix']
    tshape['CurveRepresentation']=curve
    del tshape['this']
    del tshape['ShapeType']
    del tshape['className']
    del tshape['Free']
    del tshape['Locked']
    del tshape['Modified']
    del tshape['Infinite']
    del tshape['Convex']
    del tshape['Checked']
    del tshape['Closed']
    del tshape['NbChildren']
    del curve['Location']
    return tshape
def convert_wire(w):
    obj=json.loads('{'+str(w)+'}')
    tshape=obj['TShape']
    location=obj['Location']
    del tshape['this']
    del tshape['ShapeType']
    del tshape['className']
    del tshape['Free']
    del tshape['Locked']
    del tshape['Modified']
    del tshape['Infinite']
    obj['TShape']=tshape
    return obj

def convert_face(f):
    obj=json.loads('{'+str(f)+'}')
    tshape=obj["TShape"]
    obj['Matrix']=obj['Location']['Transformation']['Matrix']
    obj['Surface']=tshape['Surface']
    del obj['Surface']['className']
    del obj['className']
    del tshape['this']
    del tshape['ShapeType']
    del tshape['Checked']
    del tshape['Closed']
    del tshape['NbChildren']
    del tshape['Location']
    del obj["TShape"]
    del obj['Location']
    return str(obj)

def convert_solid(s):
    obj=json.loads('{'+str(s)+'}')
    del obj['className']
    tshape=obj["TShape"]
    obj['NumberChildren'] = tshape['NbChildren']
    obj['Matrix']=obj['Location']['Transformation']['Matrix']
    del obj["TShape"]
    del obj['Location']
    if obj['NumberChildren'] != 1 :
        obj['children'] = [convert_solid(i) for i in s.solids]
    else:
        obj['faces'] = [convert_face(f) for f in s.faces]
        obj['vertices'] = [convert_vertex(f) for f in s.vertices]
        obj['edges'] = [convert_edge(e) for e in s.edges]
        obj['wires'] = [convert_wire(w) for w in s.wires]
    return obj

def SolidToJson(s):
    return json.dumps(convert_solid(s))

assert([convert_solid(all)==convert_solid(i) for i in all.solids]==[False]*len(all.solids))
assert(SolidToJson(all) != SolidToJson(b))
