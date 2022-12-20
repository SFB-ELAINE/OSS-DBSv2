from src.fastfouriertransform import FFTWave, Signal, FastFourierTransform
import numpy as np


class MockSignal(Signal):
    def __init__(self) -> None:
        self.frequency = 50

    def generate_samples(self, sample_spacing: float) -> np.ndarray:
        n_samples = int(1 / (self.frequency * sample_spacing))
        pulse_width = int(0.00078 * n_samples)
        padding_length = (n_samples - pulse_width)
        return np.concatenate(([1] * pulse_width, [0] * padding_length))


class TestFFTSineWave:

    def test_greater_than_amplitude(self):
        sine1 = FFTWave(amplitude=1.0, frequency=1.0)
        sine2 = FFTWave(amplitude=2.0, frequency=1.0)
        assert sine1 < sine2 and not sine1 > sine2

    def test_greater_than_amplitude_cpmplex(self):
        sine1 = FFTWave(amplitude=1.0 + 0.0j, frequency=1.0)
        sine2 = FFTWave(amplitude=3.0 + 4.0j, frequency=1.0)
        assert sine1 < sine2 and not sine1 > sine2

    def test_greater_than_frequency(self):
        sine1 = FFTWave(amplitude=1.0, frequency=2.0)
        sine2 = FFTWave(amplitude=1.0, frequency=1.0)
        assert sine1 < sine2 and not sine1 > sine2

    def test_greater_than_identical(self):
        sine1 = FFTWave(amplitude=1.0, frequency=1.0)
        assert not sine1 < sine1 and not sine1 > sine1

    def test_sorted(self):
        sine1 = FFTWave(amplitude=1.0, frequency=5.0)
        sine2 = FFTWave(amplitude=2.0, frequency=2.0)
        sine3 = FFTWave(amplitude=1.0, frequency=3.0)
        desired = [sine1, sine3, sine2]
        assert sorted([sine1, sine2, sine3]) == desired


# class TestFFT:

#     def test_sine_waves(self):
#         signal = MockSignal()
#         t = np.arange(0, 1, 1e-3)
#         fft = FFT(signal=signal)
#         fft.SPACING_FACTOR = 1e3
#       #  waves = fft.octave_bands()
#         waves = fft.sine_waves()
#         print(len(waves))
#         sine = np.zeros(len(t))
#         total_A = 0
#         for wave in waves[1:250]:
#             frequency = wave.frequency
#             A = wave.amplitude
#             total_A += A
#             phase = np.angle(A)
#             sine += abs(A) * np.cos(2*np.pi*frequency * t + phase)
 

#         sine += abs(waves[0].amplitude) / 2
#         total_A += waves[0].amplitude / 2
#         print(total_A, abs(total_A), waves[0].amplitude)
#         import matplotlib.pyplot as plt
#         plt.plot(t, sine)#, t, signal.generate_samples(1e-3/50))
#         plt.show()

