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
        return {'Flags': edge_data['TShape']['Flags'],
                'Orientable': edge_data['TShape']['Orientable'],
                'Tolerance': edge_data['TShape']['Tolerance'],
                'Curve': curve_data}

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
