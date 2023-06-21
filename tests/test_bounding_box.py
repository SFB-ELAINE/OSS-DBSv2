from ossdbs.model_geometry.bounding_box import BoundingBox


class TestBoundingBox:

    def test_shape(self):
        bbox = BoundingBox(start=(0.5, 1, -2), end=(1.4, 3, 2.05))
        assert bbox.shape == (1, 2, 4)
