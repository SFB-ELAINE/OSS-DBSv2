# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from ossdbs.utils.collapse_vta import get_collapsed_VTA


class TestCollapseVTA:
    """Tests for VTA collapse functionality."""

    def test_basic_collapse(self):
        """Test basic VTA collapse with simple inputs."""
        # Create field points: 4 points around a central axis
        # Format: [x, y, z, field_val1, field_val2, field_val3, field_val4]
        field_on_points = np.array(
            [
                [2.0, 0.0, 0.0, 1.0, 0.5, 0.3, 0.1],  # Point on +x axis
                [-2.0, 0.0, 0.0, 1.0, 0.5, 0.3, 0.1],  # Point on -x axis
                [0.0, 2.0, 0.0, 1.0, 0.5, 0.3, 0.1],  # Point on +y axis
                [0.0, -2.0, 0.0, 1.0, 0.5, 0.3, 0.1],  # Point on -y axis
            ]
        )
        implantation_coordinate = np.array([0.0, 0.0, 0.0])
        lead_direction = np.array([0.0, 0.0, 1.0])  # Lead along z-axis
        lead_diam = 1.0

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        # Check that field values are preserved
        np.testing.assert_array_equal(result[:, 3:], field_on_points[:, 3:])

        # Check that points moved inward (toward the z-axis)
        # Original distance from z-axis was 2.0, should now be 2.0 - 0.5 = 1.5
        for i in range(4):
            original_dist = np.linalg.norm(field_on_points[i, :2])
            new_dist = np.linalg.norm(result[i, :2])
            expected_dist = original_dist - lead_diam / 2.0
            np.testing.assert_almost_equal(new_dist, expected_dist, decimal=10)

    def test_collapse_preserves_z_coordinate(self):
        """Test that z-coordinates are preserved when lead is along z-axis."""
        field_on_points = np.array(
            [
                [1.0, 0.0, 5.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 10.0, 1.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, -3.0, 1.0, 0.0, 0.0, 0.0],
            ]
        )
        implantation_coordinate = np.array([0.0, 0.0, 0.0])
        lead_direction = np.array([0.0, 0.0, 1.0])
        lead_diam = 0.5

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        # Z-coordinates should be unchanged
        np.testing.assert_array_almost_equal(result[:, 2], field_on_points[:, 2])

    def test_collapse_with_offset_implantation(self):
        """Test collapse with non-origin implantation coordinate."""
        field_on_points = np.array(
            [
                [12.0, 10.0, 5.0, 1.0, 2.0, 3.0, 4.0],
            ]
        )
        implantation_coordinate = np.array([10.0, 10.0, 5.0])
        lead_direction = np.array([0.0, 0.0, 1.0])
        lead_diam = 1.0

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        # Point should move 0.5 units toward the lead axis
        # Original: (12, 10, 5), Lead axis at x=10, y=10
        # Distance from axis: 2.0, should become 1.5
        expected_x = 10.0 + 1.5  # Moved inward from 12 toward 10
        np.testing.assert_almost_equal(result[0, 0], expected_x, decimal=10)
        np.testing.assert_almost_equal(result[0, 1], 10.0, decimal=10)  # y unchanged

    def test_collapse_with_angled_lead(self):
        """Test collapse with lead not along a principal axis."""
        field_on_points = np.array(
            [
                [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            ]
        )
        implantation_coordinate = np.array([0.0, 0.0, 0.0])
        # Lead at 45 degrees in x-z plane
        lead_direction = np.array([1.0, 0.0, 1.0])
        lead_diam = 0.5

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        # Result should be valid (not NaN or Inf)
        assert np.all(np.isfinite(result))
        # Field values should be preserved
        np.testing.assert_array_equal(result[:, 3:], field_on_points[:, 3:])

    def test_collapse_output_shape(self):
        """Test that output shape matches input shape."""
        n_points = 100
        field_on_points = np.random.rand(n_points, 7)
        field_on_points[:, :3] += 5  # Offset so points are away from origin
        implantation_coordinate = np.array([0.0, 0.0, 0.0])
        lead_direction = np.array([0.0, 0.0, 1.0])
        lead_diam = 1.0

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        assert result.shape == field_on_points.shape

    def test_collapse_with_zero_diameter(self):
        """Test that zero diameter results in unchanged coordinates."""
        field_on_points = np.array(
            [
                [2.0, 0.0, 0.0, 1.0, 0.5, 0.3, 0.1],
                [0.0, 2.0, 0.0, 1.0, 0.5, 0.3, 0.1],
            ]
        )
        implantation_coordinate = np.array([0.0, 0.0, 0.0])
        lead_direction = np.array([0.0, 0.0, 1.0])
        lead_diam = 0.0

        result = get_collapsed_VTA(
            field_on_points, implantation_coordinate, lead_direction, lead_diam
        )

        # With zero diameter, coordinates should be unchanged
        np.testing.assert_array_almost_equal(result[:, :3], field_on_points[:, :3])
