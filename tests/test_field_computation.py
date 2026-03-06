# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from ossdbs.utils.field_computation import (
    compute_field_magnitude,
    compute_field_magnitude_from_components,
)


class TestComputeFieldMagnitude:
    """Tests for compute_field_magnitude function."""

    def test_unit_vectors(self):
        """Test magnitude computation for unit vectors."""
        fields = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        result = compute_field_magnitude(fields)
        expected = np.array([[1.0], [1.0], [1.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_zero_vector(self):
        """Test magnitude of zero vector."""
        fields = np.array([[0.0, 0.0, 0.0]])
        result = compute_field_magnitude(fields)
        expected = np.array([[0.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_known_magnitude(self):
        """Test with vectors of known magnitude."""
        # Vector (3, 4, 0) has magnitude 5
        fields = np.array([[3.0, 4.0, 0.0]])
        result = compute_field_magnitude(fields)
        expected = np.array([[5.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_3d_vector(self):
        """Test with full 3D vector."""
        # Vector (1, 2, 2) has magnitude 3
        fields = np.array([[1.0, 2.0, 2.0]])
        result = compute_field_magnitude(fields)
        expected = np.array([[3.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_multiple_vectors(self):
        """Test with multiple vectors."""
        fields = np.array(
            [
                [1.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [1.0, 2.0, 2.0],
            ]
        )
        result = compute_field_magnitude(fields)
        expected = np.array([[1.0], [5.0], [3.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_output_shape(self):
        """Test that output shape is (N, 1)."""
        n_points = 50
        fields = np.random.rand(n_points, 3)
        result = compute_field_magnitude(fields)
        assert result.shape == (n_points, 1)

    def test_invalid_shape_1d(self):
        """Test that 1D array raises ValueError."""
        fields = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError):
            compute_field_magnitude(fields)

    def test_invalid_shape_wrong_columns(self):
        """Test that array with wrong number of columns raises ValueError."""
        fields = np.array([[1.0, 2.0]])  # Only 2 columns instead of 3
        with pytest.raises(ValueError):
            compute_field_magnitude(fields)

    def test_negative_components(self):
        """Test with negative field components."""
        fields = np.array([[-3.0, -4.0, 0.0]])
        result = compute_field_magnitude(fields)
        expected = np.array([[5.0]])  # Magnitude is always positive
        np.testing.assert_array_almost_equal(result, expected)


class TestComputeFieldMagnitudeFromComponents:
    """Tests for compute_field_magnitude_from_components function."""

    def test_real_components(self):
        """Test with real-valued components."""
        Ex = np.array([3.0])
        Ey = np.array([4.0])
        Ez = np.array([0.0])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([5.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_complex_components(self):
        """Test with complex-valued components."""
        # |1+1j|^2 = 2, so sqrt(2 + 2 + 2) = sqrt(6)
        Ex = np.array([1.0 + 1.0j])
        Ey = np.array([1.0 + 1.0j])
        Ez = np.array([1.0 + 1.0j])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([np.sqrt(6.0)])
        np.testing.assert_array_almost_equal(result, expected)

    def test_zero_components(self):
        """Test with zero components."""
        Ex = np.array([0.0])
        Ey = np.array([0.0])
        Ez = np.array([0.0])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([0.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_multiple_points(self):
        """Test with multiple points."""
        Ex = np.array([1.0, 3.0, 0.0])
        Ey = np.array([0.0, 4.0, 0.0])
        Ez = np.array([0.0, 0.0, 5.0])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([1.0, 5.0, 5.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_pure_imaginary(self):
        """Test with pure imaginary components."""
        Ex = np.array([3.0j])
        Ey = np.array([4.0j])
        Ez = np.array([0.0j])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([5.0])  # |3j|^2 + |4j|^2 = 9 + 16 = 25, sqrt = 5
        np.testing.assert_array_almost_equal(result, expected)

    def test_mixed_complex(self):
        """Test with mixed real and imaginary components."""
        # Ex = 1, Ey = 1j, Ez = 0
        # |1|^2 = 1, |1j|^2 = 1, |0|^2 = 0
        # sqrt(1 + 1 + 0) = sqrt(2)
        Ex = np.array([1.0])
        Ey = np.array([1.0j])
        Ez = np.array([0.0])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        expected = np.array([np.sqrt(2.0)])
        np.testing.assert_array_almost_equal(result, expected)

    def test_output_shape(self):
        """Test that output shape matches input shape."""
        n_points = 100
        Ex = np.random.rand(n_points) + 1j * np.random.rand(n_points)
        Ey = np.random.rand(n_points) + 1j * np.random.rand(n_points)
        Ez = np.random.rand(n_points) + 1j * np.random.rand(n_points)
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        assert result.shape == (n_points,)

    def test_result_is_real(self):
        """Test that result is always real-valued."""
        Ex = np.array([1.0 + 2.0j, 3.0 - 1.0j])
        Ey = np.array([2.0 - 1.0j, 1.0 + 4.0j])
        Ez = np.array([0.5 + 0.5j, 2.0 + 2.0j])
        result = compute_field_magnitude_from_components(Ex, Ey, Ez)
        assert np.isrealobj(result) or np.all(np.isreal(result))
