import json
import netgen


class GeometryConverter:

    def __init__(self, solid: netgen.libngpy._NgOCC.TopoDS_Shape) -> None:
        self.__solid = solid

    def to_json(self,
                path: str) -> None:
        with open(path, "w") as file:
            file.write(json.dumps(self.to_dictionary()))

    def to_dictionary(self) -> dict:
        solid_data = self.__solid_to_dict(self.__solid)
        edges_data = [self.__edge_to_dict(e) for e in self.__solid.edges]
        faces_data = [self.__face_to_dict(face) for face in self.__solid.faces]
        vertices_data = [self.__vertex_to_dict(v)
                         for v in self.__solid.vertices]
        wires_data = [self.__wire_to_dict(w) for w in self.__solid.wires]

        if solid_data['NbChildren'] == 1:
            return {'solid': solid_data,
                    'edges': edges_data,
                    'faces': faces_data,
                    'vertices': vertices_data,
                    'wires': wires_data}

        children = [GeometryConverter(solid).to_dictionary()
                    for solid in self.__solid.solids]
        return {'solid': solid_data, 'children': children}

    @staticmethod
    def __edge_to_dict(edge: netgen.libngpy._NgOCC.Edge) -> dict:
        edge_data = json.loads("{" + str(edge) + "}")
        curve_data = edge_data['TShape']['CurveRepresentation']
        uv1 = curve_data['UV1'] if 'UV1' in  curve_data.keys() else None
        uv2 = curve_data['UV2'] if 'UV2' in curve_data.keys() else None
        pcurve = curve_data['PCurve'] if  'PCurve' in curve_data.keys() else None
        surface = curve_data['Surface'] if 'Surface' in curve_data.keys() else None
        last = curve_data['Last'] if 'Last' in curve_data.keys() else None
        first = curve_data['First'] if 'First' in curve_data.keys() else None
        ccurve = curve_data['curve'] if 'curve' in curve_data.keys() else None
        curve = {'className': curve_data['className'],
                 'Matrix': curve_data['Location']['Transformation']['Matrix']}
        if uv1 != None: curve.setdefault('UV1', uv1)
        if uv2 != None: curve.setdefault('UV2', uv2)
        if surface != None: curve.setdefault('Surface', surface)
        if pcurve != None: curve.setdefault('PCurve', pcurve)
        if last != None: curve.setdefault('Last', last)
        if first != None: curve.setdefault('First', last)
        if ccurve != None: curve.setdefault('curve', ccurve)
        return {'Flags': edge_data['TShape']['Flags'],
                'Orientable': edge_data['TShape']['Orientable'],
                'Tolerance': edge_data['TShape']['Tolerance'],
                'Curve': curve}

    @staticmethod
    def __face_to_dict(face: netgen.libngpy._NgOCC.Face) -> dict:
        face_data = json.loads('{' + str(face) + '}')
        return {'Orient': face_data['Orient'],
                'Matrix': face_data['Location']['Transformation']['Matrix'],
                'Surface': face_data["TShape"]['Surface']}

    @staticmethod
    def __solid_to_dict(solid: netgen.libngpy._NgOCC.TopoDS_Shape) -> dict:
        solid_data = json.loads('{' + str(solid) + '}')
        return {'Orient': solid_data['Orient'],
                'Matrix': solid_data['Location']['Transformation']['Matrix'],
                'NbChildren': solid_data['TShape']['NbChildren']}

    @staticmethod
    def __vertex_to_dict(vertex: netgen.libngpy._NgOCC.Vertex) -> dict:
        vertex_data = json.loads('{' + str(vertex) + '}')
        return {'transformation': vertex_data['Location']['Transformation'],
                'coordinates': vertex_data['TShape']['Pnt']['gp_Pnt']}

    @staticmethod
    def __wire_to_dict(wire: netgen.libngpy._NgOCC.Wire) -> dict:
        wire_data = json.loads('{' + str(wire) + '}')
        return {'className': wire_data['className'],
                'ShapeType': wire_data['TShape']['ShapeType'],
                'NbChildren': wire_data['TShape']['NbChildren'],
                'Flags': wire_data['TShape']['Flags'],
                'Location': wire_data['Location'],
                'Orient': wire_data['Orient']}


from netgen.occ import Sphere, Pnt, Vec


def test_foo():
    sp1 = Sphere(Pnt(0, 0, 0), 0.2)
    sp2 = Sphere(Pnt(0.5, 0, 0), 0.2)
    sp3 = Sphere(Pnt(0, 0, 0.5), 0.2)
    sp4 = Sphere(Pnt(0, 0.2, 0.7), 0.2)

    vector = Vec(1, 1, 1)
    all = sp1 + sp2 + sp3 + sp4
    compound = all.Move(vector)
    result = GeometryConverter().to_dictionary(compound)
    with open("sample.json", "w") as outfile:
        outfile.write(json.dumps(result))
