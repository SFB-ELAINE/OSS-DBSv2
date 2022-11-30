from src.signals import RectangleSignal, TriangleSignal, TrapzoidSignal
import numpy as np


class TestRectangleSignal:

    def test_generate_samples_0Hz(self):
        signal = RectangleSignal(frequency=0, pulse_width=0.5)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_1Hz(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.5)
        desired = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_3Hz(self):
        signal = RectangleSignal(frequency=3, pulse_width=0.5)
        desired = np.array([1, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_no_spacing(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_minimum_spacing(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_to_large_spacing(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=2)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_spacing(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.2)
        desired = np.append([1] * 20, [0] * 80)
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = RectangleSignal(frequency=1, pulse_width=0.0)
        desired = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_max_pulse_width(self):
        signal = RectangleSignal(frequency=1, pulse_width=1.0)
        desired = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)


class TestTriangleSignal:

    def test_generate_samples_0Hz(self):
        signal = TriangleSignal(frequency=0, pulse_width=0.5)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_1Hz(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.5)
        desired = np.array([1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_3Hz(self):
        signal = TriangleSignal(frequency=3, pulse_width=0.5)
        desired = np.array([1, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_no_spacing(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_minimum_spacing(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=1.0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_to_large_spacing(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.2)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=2.0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_spacing(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.2)
        ramp = np.arange(0.1, 1, 0.1)
        desired = np.concatenate((ramp, [1], np.flip(ramp), [0] * 81))
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.0)
        desired = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_max_pulse_width(self):
        signal = TriangleSignal(frequency=1, pulse_width=2.0)
        desired = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)


class TestTrapezoidSignal:

    def test_generate_samples_0Hz(self):
        signal = TrapzoidSignal(frequency=0, pulse_width=0.5, top_width=0.3)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_1Hz(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.5, top_width=0.3)
        desired = np.array([0.5, 1, 1, 1, 0.5, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_3Hz(self):
        signal = TrapzoidSignal(frequency=3, pulse_width=0.5, top_width=0.3)
        desired = np.array([1, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_no_spacing(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.5, top_width=0.3)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_to_large_spacing(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.5, top_width=0.3)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=2.0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_minimum_spacing(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.2, top_width=0.3)
        desired = np.array([0])
        actual = signal.generate_samples(sample_spacing=0)
        np.testing.assert_allclose(actual, desired)

    def test_generate_samples_spacing(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.2, top_width=0.1)
        ramp = np.arange(1/6, 1, 1/6)
        desired = np.concatenate((ramp, [1] * 8, np.flip(ramp), [0] * 80))
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, desired)

    def test_generate_no_pulse_width(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.0, top_width=0.0)
        desired = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)

    def test_generate_max_pulse_width(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=2.0, top_width=0.1)
        desired = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 0.8, 0.6, 0.4, 0.2])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, desired)
