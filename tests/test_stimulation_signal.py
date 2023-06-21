from ossdbs.stimulation_signals.rectangle_signal import RectangleSignal
from ossdbs.stimulation_signals.triangle_signal import TriangleSignal
from ossdbs.stimulation_signals.trapezoid_signal import TrapezoidSignal
import numpy as np
import pytest


class TestRectangleStimulationSignal:

    def test_generate_samples(self):
        signal = RectangleSignal(frequency=1,
                                 pulse_width=0.3,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        signal = RectangleSignal(frequency=1,
                                 pulse_width=0.3,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = RectangleSignal(frequency=1,
                                 pulse_width=0.0,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = RectangleSignal(frequency=1,
                                 pulse_width=0.3,
                                 counter_pulse_width=0.0,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            signal = RectangleSignal(frequency=0,
                                     pulse_width=0.3,
                                     counter_pulse_width=0.6,
                                     inter_pulse_width=0.1)

    def test_generate_over_length(self):
        signal = RectangleSignal(frequency=1,
                                 pulse_width=0.3,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.3)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, 0, 0, -0.5, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)


class TestTriangleStimulationSignal:

    def test_generate_samples(self):
        signal = TriangleSignal(frequency=1,
                                pulse_width=0.3,
                                counter_pulse_width=0.6,
                                inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, -1/6, -1/3, -0.5, -1/3, -1/6, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        signal = TriangleSignal(frequency=1,
                                pulse_width=0.3,
                                counter_pulse_width=0.6,
                                inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = TriangleSignal(frequency=1,
                                pulse_width=0.0,
                                counter_pulse_width=0.6,
                                inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = TriangleSignal(frequency=1,
                                pulse_width=0.3,
                                counter_pulse_width=0.0,
                                inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            signal = TriangleSignal(frequency=0,
                                    pulse_width=0.3,
                                    counter_pulse_width=0.6,
                                    inter_pulse_width=0.1)

    def test_generate_over_length(self):
        signal = TriangleSignal(frequency=1,
                                pulse_width=0.3,
                                counter_pulse_width=0.6,
                                inter_pulse_width=0.3)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, -1/6, -1/3, -0.5, -1/3]
        np.testing.assert_equal(actual, desired)


class TestTrapezoidStimulationSignal:

    def test_generate_samples(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.3,
                                 top_width=0.1,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, -1 / 6, -1 / 3, -0.5, -0.5, -1 / 3, -1 / 6]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.3,
                                 top_width=0.1,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    """
    # TODO fix
    def test_generate_no_pulse_width(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.0,
                                 top_width=0.1,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)
    """

    def test_generate_no_pulse_top_width(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.3,
                                 top_width=0.0,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.3,
                                 top_width=0.1,
                                 counter_pulse_width=0.0,
                                 inter_pulse_width=0.1)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        with pytest.raises(ValueError):
            signal = TrapezoidSignal(frequency=0,
                                     pulse_width=0.3,
                                     top_width=0.1,
                                     counter_pulse_width=0.6,
                                     inter_pulse_width=0.1)

    def test_generate_over_length(self):
        signal = TrapezoidSignal(frequency=1,
                                 pulse_width=0.3,
                                 top_width=0.1,
                                 counter_pulse_width=0.6,
                                 inter_pulse_width=0.3)
        actual = signal._generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, -1 / 6, -1 / 3, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)
