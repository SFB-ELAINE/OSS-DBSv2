from typing import ClassVar

import numpy as np
import pytest

from ossdbs.stimulation_signals.rectangle_signal import RectangleSignal
from ossdbs.stimulation_signals.trapezoid_signal import TrapezoidSignal
from ossdbs.stimulation_signals.triangle_signal import TriangleSignal
from ossdbs.stimulation_signals.utilities import (
    adjust_cutoff_frequency,
    get_indices_in_octave_band,
    get_maximum_octave_band_index,
    get_minimum_octave_band_index,
    get_octave_band_indices,
    get_timesteps,
    reconstruct_time_signals,
)


class TestFFT:
    TESTDATA: ClassVar[list[tuple[str, float, float, float, float, float]]] = [
        # signal_type, frequency, pulse_width, counter_pulse_width,
        # inter_pulse_width, cutoff_frequency
        ("Rectangle", 130.0, 60e-6, 120e-6, 120e-6, 1e5),
        ("Triangle", 130.0, 60e-6, 120e-6, 120e-6, 1e5),
        ("Trapezoid", 130.0, 60e-6, 120e-6, 120e-6, 1e5),
    ]

    @pytest.mark.parametrize(
        "signal_type, frequency, pulse_width, counter_pulse_width,\
        inter_pulse_width, cutoff_frequency",
        TESTDATA,
    )
    def test_time_domain_signal(
        self,
        signal_type,
        frequency,
        pulse_width,
        counter_pulse_width,
        inter_pulse_width,
        cutoff_frequency,
    ):
        adj_cutoff_frequency = adjust_cutoff_frequency(
            2.0 * cutoff_frequency, frequency
        )
        dt = 1.0 / adj_cutoff_frequency
        timesteps = int(adj_cutoff_frequency / frequency)

        if signal_type == "Rectangle":
            signal = RectangleSignal(
                frequency, pulse_width, inter_pulse_width, counter_pulse_width
            )
        elif signal_type == "Triangle":
            signal = TriangleSignal(
                frequency, pulse_width, inter_pulse_width, counter_pulse_width
            )
        elif signal_type == "Trapezoid":
            top_width = 30e-6
            signal = TrapezoidSignal(
                frequency,
                pulse_width,
                inter_pulse_width,
                top_width,
                counter_pulse_width,
            )

        original_signal = signal.get_time_domain_signal(dt, timesteps)
        _, fft_signal, signal_length = signal.get_fft_spectrum(cutoff_frequency)
        _, retrieved_signal = signal.retrieve_time_domain_signal(
            fft_signal, cutoff_frequency, signal_length
        )
        tolerance = 1e-5
        np.testing.assert_allclose(original_signal, retrieved_signal, atol=tolerance)


class TestRectangleStimulationSignal:
    def test_generate_samples(self):
        signal = RectangleSignal(
            frequency=0.5,
            pulse_width=0.3,
            counter_pulse_width=0.5,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 1, 1, 1, 0, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = RectangleSignal(
                frequency=0.5,
                pulse_width=0.3,
                counter_pulse_width=0.5,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=15)

    def test_generate_no_pulse_width(self):
        signal = RectangleSignal(
            frequency=0.5,
            pulse_width=0.0,
            counter_pulse_width=0.5,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = RectangleSignal(
            frequency=0.5,
            pulse_width=0.3,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            RectangleSignal(
                frequency=0,
                pulse_width=0.3,
                counter_pulse_width=0.5,
                inter_pulse_width=0.1,
            )

    def test_generate_over_length(self):
        signal = RectangleSignal(
            frequency=0.5,
            pulse_width=0.3,
            counter_pulse_width=0.5,
            inter_pulse_width=0.3,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 1, 1, 1, 0, 0, 0, -0.5, -0.5, -0.5, -0.5, -0.5, 0]
        np.testing.assert_equal(actual, desired)


class TestTriangleStimulationSignal:
    def test_generate_samples(self):
        signal = TriangleSignal(
            frequency=0.5,
            pulse_width=0.2,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 1, 0, 0, -1 / 6, -1 / 3, -0.5, -1 / 3, -1 / 6, 0, 0, 0, 0]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = TriangleSignal(
                frequency=0.5,
                pulse_width=0.5,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=15)

    def test_generate_no_pulse_width(self):
        signal = TriangleSignal(
            frequency=0.5,
            pulse_width=0.0,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, -1 / 6, -1 / 3, -0.5, -1 / 3, -1 / 6, 0, 0, 0, 0, 0, 0, 0, 0]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_generate_no_counter_pulse_width(self):
        signal = TriangleSignal(
            frequency=0.5,
            pulse_width=0.3,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 0, 2 / 3, 2 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            TriangleSignal(
                frequency=0,
                pulse_width=0.3,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
            )

    def test_generate_over_length(self):
        signal = TriangleSignal(
            frequency=0.5,
            pulse_width=0.2,
            counter_pulse_width=0.6,
            inter_pulse_width=0.2,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 1, 0, 0, 0, -1 / 6, -1 / 3, -0.5, -1 / 3, -1 / 6, 0, 0, 0]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)


class TestTrapezoidStimulationSignal:
    def test_generate_samples(self):
        signal = TrapezoidSignal(
            frequency=0.5,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, -0.25, -0.5, -0.5, -0.25, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = TrapezoidSignal(
                frequency=0.5,
                pulse_width=0.3,
                top_width=0.1,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=15)

    def test_generate_no_pulse_width(self):
        signal = TrapezoidSignal(
            frequency=0.5,
            pulse_width=0.0,
            top_width=0.1,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, -0.25, -0.5, -0.5, -0.25, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_generate_no_counter_pulse_width(self):
        signal = TrapezoidSignal(
            frequency=0.5,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            TrapezoidSignal(
                frequency=0,
                pulse_width=0.3,
                top_width=0.1,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )

    def test_generate_over_length(self):
        signal = TrapezoidSignal(
            frequency=0.5,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.6,
            inter_pulse_width=0.3,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=15)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -0.25, -0.5, -0.5, -0.25, 0]
        np.testing.assert_equal(actual, desired)


class TestAdjustCutoffFrequency:
    """Tests for adjust_cutoff_frequency function."""

    def test_already_multiple(self):
        """Test when cutoff is already a multiple of frequency."""
        result = adjust_cutoff_frequency(1300.0, 130.0)
        assert result == 1300.0

    def test_needs_adjustment(self):
        """Test when cutoff needs adjustment."""
        result = adjust_cutoff_frequency(1350.0, 130.0)
        assert result == 1300.0

    def test_small_remainder(self):
        """Test with small remainder."""
        result = adjust_cutoff_frequency(1301.0, 130.0)
        assert result == 1300.0

    def test_large_frequency(self):
        """Test with large cutoff frequency."""
        result = adjust_cutoff_frequency(100000.0, 130.0)
        expected = 100000.0 - (100000.0 % 130.0)
        assert result == expected


class TestGetTimesteps:
    """Tests for get_timesteps function."""

    def test_basic_timesteps(self):
        """Test basic timestep generation."""
        cutoff_frequency = 1000.0
        base_frequency = 100.0
        n_frequencies = 10
        result = get_timesteps(cutoff_frequency, base_frequency, n_frequencies)
        assert len(result) == n_frequencies
        assert result[0] == 0.0

    def test_timestep_spacing(self):
        """Test that timesteps are evenly spaced."""
        cutoff_frequency = 1000.0
        base_frequency = 100.0
        n_frequencies = 20
        result = get_timesteps(cutoff_frequency, base_frequency, n_frequencies)
        diffs = np.diff(result)
        np.testing.assert_array_almost_equal(diffs, diffs[0] * np.ones_like(diffs))


class TestReconstructTimeSignals:
    """Tests for reconstruct_time_signals function."""

    def test_simple_reconstruction(self):
        """Test reconstruction of simple frequency domain signal."""
        # Create a simple frequency domain signal
        signal_length = 16
        freq_domain_signal = np.zeros(signal_length // 2 + 1, dtype=complex)
        freq_domain_signal[0] = signal_length  # DC component
        result = reconstruct_time_signals(freq_domain_signal, signal_length)
        # DC signal should result in constant time signal
        np.testing.assert_array_almost_equal(result, np.ones(signal_length))

    def test_reconstruction_preserves_length(self):
        """Test that reconstruction preserves signal length."""
        signal_length = 128
        freq_domain_signal = np.random.rand(
            signal_length // 2 + 1
        ) + 1j * np.random.rand(signal_length // 2 + 1)
        result = reconstruct_time_signals(freq_domain_signal, signal_length)
        assert len(result) == signal_length


class TestOctaveBandFunctions:
    """Tests for octave band related functions."""

    def test_get_octave_band_indices(self):
        """Test octave band index generation."""
        frequencies = np.arange(0, 17)  # 0 to 16
        result = get_octave_band_indices(frequencies)
        expected = np.array([1, 2, 4, 8, 16])
        np.testing.assert_array_equal(result, expected)

    def test_get_minimum_octave_band_index(self):
        """Test minimum octave band index calculation."""
        result = get_minimum_octave_band_index(8)
        expected = int(np.round(8 / np.sqrt(2)))
        assert result == expected

    def test_get_maximum_octave_band_index(self):
        """Test maximum octave band index calculation."""
        result = get_maximum_octave_band_index(8)
        expected = int(np.round(8 * np.sqrt(2)))
        assert result == expected

    def test_get_indices_in_octave_band_basic(self):
        """Test basic octave band indices retrieval."""
        freq_idx = 4
        frequency_indices = np.array([1, 2, 4, 8])
        cutoff_frequency_index = 16
        result = get_indices_in_octave_band(
            freq_idx, frequency_indices, cutoff_frequency_index
        )
        assert isinstance(result, list | np.ndarray)
        assert len(result) > 0

    def test_get_indices_in_octave_band_first_band(self):
        """Test octave band indices for first frequency."""
        freq_idx = 1
        frequency_indices = np.array([1, 2, 4, 8])
        cutoff_frequency_index = 16
        result = get_indices_in_octave_band(
            freq_idx, frequency_indices, cutoff_frequency_index
        )
        assert freq_idx in result or len(result) > 0

    def test_get_indices_in_octave_band_invalid_frequency(self):
        """Test that invalid frequency raises error."""
        freq_idx = 5  # Not in the list
        frequency_indices = np.array([1, 2, 4, 8])
        cutoff_frequency_index = 16
        with pytest.raises(ValueError):
            get_indices_in_octave_band(
                freq_idx, frequency_indices, cutoff_frequency_index
            )
