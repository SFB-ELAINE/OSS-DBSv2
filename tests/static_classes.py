from json import JSONEncoder 
import dataclasses

# Dataclass JSONEncoder, special behaviour if object is a dataclass object
class EnhancedJSONEncoder(JSONEncoder):
        def default(self, o):
            if dataclasses.is_dataclass(o):
                return dataclasses.asdict(o)
            return super().default(o)