import netgen
from static_classes import EnhancedJSONEncoder
import json
from create_json import convert_solid


class ElectrodeTestBase():
    def electrode_to_dict(self, geometry:netgen.libngpy._NgOCC.TopoDS_Shape):
        # Encode dataclasses as string and decode it as dict
        return json.loads(json.dumps(convert_solid(geometry),cls=EnhancedJSONEncoder))
        
