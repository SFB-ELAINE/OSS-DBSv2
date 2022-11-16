from src.signals import RectangleSignal, TriangleSignal, TrapzoidSignal
import numpy as np


class TestRectangleSignal:

    def test_generate_samples_1Hz(self):
        signal = RectangleSignal(frequency=1, pulswidth=0.5)
        expected = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
        assert np.all(signal.generate_samples(sample_spacing=0.1) == expected)

    def test_generate_samples_3Hz(self):
        signal = RectangleSignal(frequency=3, pulswidth=0.5)
        expected = np.array([1, 0, 0])
        assert np.all(signal.generate_samples(sample_spacing=0.1) == expected)

    def test_generate_samples_spacing(self):
        signal = RectangleSignal(frequency=1, pulswidth=0.2)
        expected = np.append([1] * 20, [0] * 80)
        assert np.all(signal.generate_samples(sample_spacing=0.01) == expected)


class TestTriangleSignal:

    def test_generate_samples_1Hz(self):
        signal = TriangleSignal(frequency=1, pulswidth=0.5)
        expected = np.array([1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, expected)

    def test_generate_samples_3Hz(self):
        signal = TriangleSignal(frequency=3, pulswidth=0.5)
        expected = np.array([1, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, expected)

    def test_generate_samples_spacing(self):
        signal = TriangleSignal(frequency=1, pulswidth=0.2)
        right_ramp = np.arange(1, 0, -0.1)
        left_ramp = np.flip(right_ramp)
        expected = np.concatenate((left_ramp[:-1], right_ramp, [0] * 81))
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, expected)


class TestTrapezoidSignal:

    def test_generate_samples_1Hz(self):
        signal = TrapzoidSignal(frequency=1, pulswidth=0.5, top_width=0.3)
        expected = np.array([0.5, 1, 1, 1, 0.5, 0, 0, 0, 0, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, expected)

    def test_generate_samples_3Hz(self):
        signal = TrapzoidSignal(frequency=3, pulswidth=0.5, top_width=0.3)
        expected = np.array([1, 1, 0])
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, expected)

    # def test_generate_samples_spacing(self):
    #     signal = TrapzoidSignal(frequency=1, pulswidth=0.2)
    #     right_ramp = np.arange(1, 0, -0.1)
    #     left_ramp = np.flip(right_ramp)
    #     expected = np.concatenate((left_ramp[:-1], right_ramp, [0] * 81))
    #     actual = signal.generate_samples(sample_spacing=0.01)
    #     np.testing.assert_allclose(actual, expected)
