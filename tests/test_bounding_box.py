import pytest

from ossdbs.model_geometry.bounding_box import BoundingBox


class TestBoundingBox:
    @pytest.fixture
    def bbox(self):
        return BoundingBox(start=(0.5, 1, -2), end=(1.4, 3, 2.05))

    def test_shape(self, bbox):
        assert bbox.shape == (1, 2, 4)

    def test_intersection(self, bbox):
        bbox2 = BoundingBox(start=(1.0, -1.0, 2.0), end=(1.0, 3.0, 5.0))
        intersect = bbox.intersection(bbox2)

        assert intersect.start == (1.0, 1.0, 2.0) and intersect.end == (1.0, 3.0, 2.05)

    def test_points(self, bbox):
        desired = [
            (-1.0, 1.0, -3.0),
            (-1.0, 1.0, -1.0),
            (-1.0, 1.0, 1.0),
            (1.0, 1.0, -3.0),
            (1.0, 1.0, -1.0),
            (1.0, 1.0, 1.0),
        ]

        points = bbox.points(offset=(1.0, 1.0, 1.0), voxel_size=(2.0, 2.0, 2.0))
        assert points == desired
