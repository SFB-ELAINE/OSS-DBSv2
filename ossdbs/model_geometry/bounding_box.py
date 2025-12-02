# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np


class BoundingBox:
    """Represents a cuboid aligned with the cartesian axis.

    Attributes
    ----------
    start : tuple
        Start point (x, y, z) of bounding box.

    end : tuple
        End point (x, y, z) of bounding box.
    """

    def __init__(self, start: tuple = (0, 0, 0), end: tuple = (0, 0, 0)) -> None:
        self.start = start
        self.end = end

    @property
    def shape(self) -> tuple[int]:
        """Return the integer distances between start and end point of
        bounding box.

        Returns
        -------
        tuple of int
        """
        return tuple(np.round(np.subtract(self.end, self.start)).astype(int))

    def intersection(self, other: "BoundingBox") -> "BoundingBox":
        """Returns the overlapping volume of this and another bounding boxes.

        Parameters
        ----------
        other: BoundingBox
           Second bounding box

        Returns
        -------
        BoundingBox
        """
        x_s = max(self.start[0], other.start[0])
        y_s = max(self.start[1], other.start[1])
        z_s = max(self.start[2], other.start[2])
        x_e = min(self.end[0], other.end[0])
        y_e = min(self.end[1], other.end[1])
        z_e = min(self.end[2], other.end[2])

        if x_s > x_e or y_s > y_e or z_s > z_e:
            return BoundingBox(start=(0, 0, 0), end=(0, 0, 0))

        return BoundingBox(start=(x_s, y_s, z_s), end=(x_e, y_e, z_e))

    def points(self, offset: tuple, voxel_size: tuple) -> list[tuple]:
        """Generate point coordinates matrix.

        Parameters
        ----------
        offset : tuple
            Shift of the point matrix along the cartesian coordinates
            (x, y, z).

        voxel_size : tuple
            Distances (x, y, z) between adjacent points along the cartesian
            coordinates.
        """
        start, end = self.start, self.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        shape = end_index - start_index
        x_s, y_s, z_s = start_index * voxel_size + offset
        x_values = x_s + np.arange(shape[0]) * voxel_size[0]
        y_values = y_s + np.arange(shape[1]) * voxel_size[1]
        z_values = z_s + np.arange(shape[2]) * voxel_size[2]
        return [(x, y, z) for x in x_values for y in y_values for z in z_values]

    def __eq__(self, other: "BoundingBox") -> bool:
        """Copy other bounding box."""
        return self.start == other.start and self.end == other.end
