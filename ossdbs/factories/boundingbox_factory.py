
from ossdbs.bounding_box import BoundingBox
import numpy as np


class BoundingBoxFactory:

    @classmethod
    def create(cls, box_parameters: dict) -> BoundingBox:
        input_s = box_parameters['Shape']
        input_c = box_parameters['Center']
        shape = (input_s['x[mm]'], input_s['y[mm]'], input_s['z[mm]'])
        center = (input_c['x[mm]'], input_c['y[mm]'], input_c['z[mm]'])
        start = center - np.divide(shape, 2)
        end = start + shape
        return BoundingBox(tuple(start), tuple(end))
