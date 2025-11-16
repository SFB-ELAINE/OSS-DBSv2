from typing import ClassVar

import numpy as np
import pytest

from ossdbs.stimulation_signals.rectangle_signal import RectangleSignal
from ossdbs.stimulation_signals.trapezoid_signal import TrapezoidSignal
from ossdbs.stimulation_signals.triangle_signal import TriangleSignal
from ossdbs.stimulation_signals.utilities import adjust_cutoff_frequency


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
