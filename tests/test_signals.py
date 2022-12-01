from src.signals import RectangleSignal, TriangleSignal, TrapzoidSignal
import numpy as np
import pytest


class TestRectangleSignal:

    @pytest.mark.parametrize('frequency, samples',
                             [(0, [0]),
                              (1, [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]),
                              (3, [1, 0, 0])])
    def test_generate_samples_frequency(self, frequency, samples):
        signal = RectangleSignal(frequency=frequency, pulse_width=0.5)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_equal(actual, samples)

    @pytest.mark.parametrize('spacing, samples',
                             [(0, [0]),
                              (1, [0]),
                              (2, [0]),
                              (0.01, np.append([1] * 50, [0] * 50))])
    def test_generate_samples_sample_spacing(self, spacing, samples):
        signal = RectangleSignal(frequency=1, pulse_width=0.5)
        actual = signal.generate_samples(sample_spacing=spacing)
        np.testing.assert_equal(actual, samples)

    @pytest.mark.parametrize('pulse_width, samples',
                             [(0, [0] * 10),
                              (1, [1] * 10),
                              (0.2, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0])])
    def test_generate_samples_pulse_width(self, pulse_width, samples):
        signal = RectangleSignal(frequency=1, pulse_width=pulse_width)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_equal(actual, samples)


class TestTriangleSignal:

    @pytest.mark.parametrize('frequency, samples',
                             [(0, [0]),
                              (1, [1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0]),
                              (3, [1, 0, 0])])
    def test_generate_samples_frequency(self, frequency, samples):
        signal = TriangleSignal(frequency=frequency, pulse_width=0.5)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_equal(actual, samples)

    @pytest.mark.parametrize('spacing, samples',
                             [(0, [0]),
                              (1, [0]),
                              (2, [0])])
    def test_generate_samples_sample_spacing_no_signal(self, spacing, samples):
        signal = RectangleSignal(frequency=1, pulse_width=0.5)
        actual = signal.generate_samples(sample_spacing=spacing)
        np.testing.assert_equal(actual, samples)

    def test_generate_samples_spacing_(self):
        signal = TriangleSignal(frequency=1, pulse_width=0.2)
        ramp = np.arange(0.1, 1, 0.1)
        desired = np.concatenate((ramp, [1], np.flip(ramp), [0] * 81))
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, desired)

    @pytest.mark.parametrize('pulse_width, samples',
                             [(0, [0] * 10),
                              (1, [.2, .4, .6, .8, 1., .8, .6, .4, .2, 0.]),
                              (0.5, [1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0])])
    def test_generate_pulse_width(self, pulse_width, samples):
        signal = TriangleSignal(frequency=1, pulse_width=pulse_width)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, samples)


class TestTrapezoidSignal:

    @pytest.mark.parametrize('frequency, samples',
                             [(0, [0]),
                              (1, [1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0]),
                              (3, [1, 0, 0])])
    def test_generate_samples_frequency(self, frequency, samples):
        signal = TrapzoidSignal(frequency=frequency,
                                pulse_width=0.5,
                                top_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_equal(actual, samples)

    @pytest.mark.parametrize('spacing, samples',
                             [(0, [0]),
                              (1, [0]),
                              (2, [0])])
    def test_generate_samples_sample_spacing_no_signal(self, spacing, samples):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.5, top_width=0.1)
        actual = signal.generate_samples(sample_spacing=spacing)
        np.testing.assert_equal(actual, samples)

    def test_generate_samples_spacing(self):
        signal = TrapzoidSignal(frequency=1, pulse_width=0.2, top_width=0.1)
        ramp = np.arange(1/6, 1, 1/6)
        desired = np.concatenate((ramp, [1] * 8, np.flip(ramp), [0] * 80))
        actual = signal.generate_samples(sample_spacing=0.01)
        np.testing.assert_allclose(actual, desired)

    @pytest.mark.parametrize('pulse_width, samples',
                             [(0, [0] * 10),
                              (1, [.2, .4, .6, .8, 1., 1., .8, .6, .4, .2]),
                              (0.5, [1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0])])
    def test_generate_samples_pulse_width(self, pulse_width, samples):
        signal = TrapzoidSignal(frequency=1,
                                pulse_width=pulse_width,
                                top_width=0.1)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, samples)

    @pytest.mark.parametrize('top_width, samples',
                             [(0, [1/3, 2/3, 1, 2/3, 1/3, 0, 0, 0, 0, 0]),
                              (1, [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]),
                              (0.25, [0.5, 1, 1, 1, 0.5, 0, 0, 0, 0, 0])])
    def test_generate_samples_top_width(self, top_width, samples):
        signal = TrapzoidSignal(frequency=1,
                                pulse_width=0.5,
                                top_width=top_width)
        actual = signal.generate_samples(sample_spacing=0.1)
        np.testing.assert_allclose(actual, samples)
