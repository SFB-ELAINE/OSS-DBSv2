
from dataclasses import dataclass


@dataclass
class Point:
    x: int
    y: int
    z: int

    def __array__(self) -> np.ndarray:
        return self.x, self.y, self.z


class Tetrahedron:
    
    def __init__(self, v1: Point, v2: Point, v3: Point, v4: Point) -> None:
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4

    def centroid(self) -> Point:
        x = (self.v1.x + self.v2.x + self.v3.x + self.v4.x) / 4.0
        y = (self.v1.y + self.v2.y + self.v3.y + self.v4.y) / 4.0
        z = (self.v1.z + self.v2.z + self.v3.z + self.v4.z) / 4.0
        return Point(x, y, z)
