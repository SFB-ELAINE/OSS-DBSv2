from typing import List,Any
from dataclasses import dataclass

@dataclass
class Vertex_data:
    Transformation: dict
    Point: tuple
    def dict(self) -> dict :
        return json.dumps({"Transformation": Transformation , "Point":Point})

@dataclass 
class CurveRepresentation:
    className : str
    First : float
    Last  : float 
    UV1 : dict  # 2D vertex
    UV2  : dict # 2D vertex
    PCurve : dict
    Surface : dict
    TransformationMatrix : List[int] # List of ints usually 

@dataclass
class Edge_data:
    Flags: int
    Orientable : int
    Tolerance: float
    CurveRepresentation : CurveRepresentation 
        
@dataclass
class Wire_data: 
    className : str
    TShape : dict
    Location : dict
    Orient :int 
    

@dataclass
class Face_data:
    Orient : int
    Matrix : List[int]
    Surface :dict

@dataclass
class Solid_data:
    Children : List[Any] # Recursive type definition
    Orient : int
    NumberChildren : int
    TransformationMatrix : List[int]
    Faces : List[Face_data]
    Vertices : List[Vertex_data]
    Edges : List[Edge_data]
    Wires : List[Wire_data]



