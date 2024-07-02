import numpy as np
import pytest

from ossdbs.stimulation_signals.rectangle_signal import RectangleSignal
from ossdbs.stimulation_signals.trapezoid_signal import TrapezoidSignal
from ossdbs.stimulation_signals.triangle_signal import TriangleSignal


class TestRectangleStimulationSignal:
    def test_generate_samples(self):
        signal = RectangleSignal(
            frequency=1,
            pulse_width=0.3,
            counter_pulse_width=0.5,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=12)
        desired = [0, 0, 0, 1, 1, 1, 0, -0.5, -0.5, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = RectangleSignal(
                frequency=1,
                pulse_width=0.3,
                counter_pulse_width=0.5,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=12)

    def test_generate_no_pulse_width(self):
        signal = RectangleSignal(
            frequency=1,
            pulse_width=0.0,
            counter_pulse_width=0.5,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=12)
        desired = [0, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = RectangleSignal(
            frequency=1,
            pulse_width=0.3,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=12)
        desired = [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
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
            frequency=1,
            pulse_width=0.3,
            counter_pulse_width=0.5,
            inter_pulse_width=0.3,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=12)
        desired = [0, 0, 0, 1, 1, 1, 0, 0, 0, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)


class TestTriangleStimulationSignal:
    def test_generate_samples(self):
        signal = TriangleSignal(
            frequency=1,
            pulse_width=0.2,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=11)
        desired = [0, 0, 0, 1, 0, 0, -1 / 6, -1 / 3, -0.5, -1 / 3, -1 / 6]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = TriangleSignal(
                frequency=1,
                pulse_width=0.3,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=11)

    def test_generate_no_pulse_width(self):
        signal = TriangleSignal(
            frequency=1,
            pulse_width=0.0,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=11)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = TriangleSignal(
            frequency=1,
            pulse_width=0.3,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=11)
        desired = [0, 0, 0, 0, 2 / 3, 2 / 3, 0, 0, 0, 0, 0]
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
            frequency=1,
            pulse_width=0.3,
            counter_pulse_width=0.6,
            inter_pulse_width=0.3,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=11)
        desired = [0, 0, 0, 0, 2 / 3, 2 / 3, 0, 0, 0, 0, -1 / 6]
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)


class TestTrapezoidStimulationSignal:
    def test_generate_samples(self):
        signal = TrapezoidSignal(
            frequency=1,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=13)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, -0.25, -0.5, -0.5, -0.25, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        with pytest.raises(ValueError):
            signal = TrapezoidSignal(
                frequency=1,
                pulse_width=0.3,
                top_width=0.1,
                counter_pulse_width=0.6,
                inter_pulse_width=0.1,
                counter_pulse_amplitude=0.5,
            )
            signal.get_time_domain_signal(dt=0, timesteps=13)

    def test_generate_no_pulse_top_width(self):
        signal = TrapezoidSignal(
            frequency=1,
            pulse_width=0.3,
            top_width=0.0,
            counter_pulse_width=0.6,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=13)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = TrapezoidSignal(
            frequency=1,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.0,
            inter_pulse_width=0.1,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=13)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]
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
            frequency=1,
            pulse_width=0.3,
            top_width=0.1,
            counter_pulse_width=0.6,
            inter_pulse_width=0.3,
            counter_pulse_amplitude=0.5,
        )
        actual = signal.get_time_domain_signal(dt=0.1, timesteps=13)
        desired = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -0.25, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)
