
from dataclasses import dataclass


@dataclass
class BoundingBox:
    start: tuple = (0., 0., 0.)
    end: tuple = (0., 0., 0.)
