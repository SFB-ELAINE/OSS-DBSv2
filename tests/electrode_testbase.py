import netgen
from tests.static_classes import EnhancedJSONEncoder
import json
from tests.create_json import convert_solid
from abc import abstractmethod


class ElectrodeTestBase():
    def electrode_to_dict(self, geometry:netgen.libngpy._NgOCC.TopoDS_Shape) -> dict:
        # Encode dataclasses as string and decode it as dict
        return json.loads(json.dumps(convert_solid(geometry),cls=EnhancedJSONEncoder))
    def load_json(self, path: str)-> dict:
        with open(path, "r") as file:
            geometry_dict=json.load(file)
        return geometry_dict
    @abstractmethod
    def test_generate_geometry_default(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        pass
        
