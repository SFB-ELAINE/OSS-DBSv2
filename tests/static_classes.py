from json import JSONEncoder 

class EnhancedJSONEncoder(JSONEncoder):
        def default(self, o):
            if dataclasses.is_dataclass(o):
                return dataclasses.asdict(o)
            return super().default(o)