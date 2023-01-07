import netgen
from static_classes import EnhancedJSONEncoder
import json
from create_json import convert_solid


class ElectrodeTestBase():
    def electrode_to_dict(self, geometry:netgen.libngpy._NgOCC.TopoDS_Shape) -> dict:
        # Encode dataclasses as string and decode it as dict
        return json.loads(json.dumps(convert_solid(geometry),cls=EnhancedJSONEncoder))
    def load_json(self, path: str)-> dict:
        with open(path, "r") as file:
            geometry_dict=json.load(file)
        return geometry_dict
        
