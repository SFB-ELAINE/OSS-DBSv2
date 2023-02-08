from ossdbs.stimulation_signal import RectangleStimulationSignal
from ossdbs.stimulation_signal import TriangleStimulationSignal
import numpy as np


class TestRectangleStimulationSignal:

    def test_generate_samples(self):
        signal = RectangleStimulationSignal(frequency=1,
                                            pulse_width=0.3,
                                            counter_pulse_width=0.6,
                                            inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        signal = RectangleStimulationSignal(frequency=1,
                                            pulse_width=0.3,
                                            counter_pulse_width=0.6,
                                            inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = RectangleStimulationSignal(frequency=1,
                                            pulse_width=0.0,
                                            counter_pulse_width=0.6,
                                            inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = RectangleStimulationSignal(frequency=1,
                                            pulse_width=0.3,
                                            counter_pulse_width=0.0,
                                            inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        signal = RectangleStimulationSignal(frequency=0,
                                            pulse_width=0.3,
                                            counter_pulse_width=0.6,
                                            inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_over_length(self):
        signal = RectangleStimulationSignal(frequency=1,
                                            pulse_width=0.3,
                                            counter_pulse_width=0.6,
                                            inter_pulse_width=0.3)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [1, 1, 1, 0, 0, 0, -0.5, -0.5, -0.5, -0.5]
        np.testing.assert_equal(actual, desired)


class TestTriangleStimulationSignal:

    def test_generate_samples(self):
        signal = TriangleStimulationSignal(frequency=1,
                                           pulse_width=0.3,
                                           counter_pulse_width=0.6,
                                           inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, -1/6, -1/3, -0.5, -1/3, -1/6, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_samples(self):
        signal = TriangleStimulationSignal(frequency=1,
                                           pulse_width=0.3,
                                           counter_pulse_width=0.6,
                                           inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = TriangleStimulationSignal(frequency=1,
                                           pulse_width=0.0,
                                           counter_pulse_width=0.6,
                                           inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_counter_pulse_width(self):
        signal = TriangleStimulationSignal(frequency=1,
                                           pulse_width=0.3,
                                           counter_pulse_width=0.0,
                                           inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, 0, 0, 0, 0]
        np.testing.assert_equal(actual, desired)

    def test_generate_no_frequency(self):
        signal = TriangleStimulationSignal(frequency=0,
                                           pulse_width=0.3,
                                           counter_pulse_width=0.6,
                                           inter_pulse_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0]
        np.testing.assert_equal(actual, desired)

    def test_generate_over_length(self):
        signal = TriangleStimulationSignal(frequency=1,
                                           pulse_width=0.3,
                                           counter_pulse_width=0.6,
                                           inter_pulse_width=0.3)
        actual = signal.generate_samples(sample_spacing=0.1)
        desired = [0.5, 1, 0.5, 0, 0, 0, -1/6, -1/3, -0.5, -1/3]
        np.testing.assert_equal(actual, desired)
