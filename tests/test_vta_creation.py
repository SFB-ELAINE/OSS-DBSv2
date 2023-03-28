# from ossdbs.vta_points import VTAPointMatrix
# import numpy as np


# class TestVTAPoints:

#     def test_coordinates_one_point(self):
#         points = VTAPointMatrix(distance=0.5,
#                                 shape=(1, 1, 1),
#                                 center=(0, 0, 0),
#                                 direction=(0, 0, 1))
#         actual = points.coordinates()
#         desired = [(0, 0, 0)]
#         np.testing.assert_equal(actual, desired)

#     def test_coordinates_two_points(self):
#         points = VTAPointMatrix(distance=0.5,
#                                 shape=(1, 2, 1),
#                                 center=(0, 0, 0),
#                                 direction=(0, 0, 1))
#         actual = points.coordinates()
#         desired = [(0, -0.25, 0), (0, 0.25, 0)]
#         np.testing.assert_equal(actual, desired)

#     def test_coordinates_center(self):
#         points = VTAPointMatrix(distance=0.5,
#                                 shape=(2, 2, 1),
#                                 center=(1, -1, 2),
#                                 direction=(0, 0, 1))
#         actual = points.coordinates()
#         desired = [(0.75, -1.25, 2.0),
#                    (0.75, -0.75, 2.0),
#                    (1.25, -1.25, 2.0),
#                    (1.25, -0.75, 2.0)]
#         np.testing.assert_equal(actual, desired)

#     def test_coordinates_distance(self):
#         points = VTAPointMatrix(distance=2,
#                                 shape=(2, 2, 1),
#                                 center=(0, 0, 0),
#                                 direction=(0, 0, 1))
#         actual = points.coordinates()
#         desired = [(-1, -1, 0),
#                    (-1,  1, 0),
#                    (1, -1, 0),
#                    (1,  1, 0)]
#         np.testing.assert_equal(actual, desired)

#     def test_coordinates_three_points(self):
#         points = VTAPointMatrix(distance=0.5,
#                                 shape=(3, 1, 1),
#                                 center=(0, 0, 0),
#                                 direction=(0, 0, 1))
#         actual = points.coordinates()
#         desired = [(-0.5, 0, 0),
#                    (0,  0, 0),
#                    (0.5, 0, 0)]
#         np.testing.assert_equal(actual, desired)
