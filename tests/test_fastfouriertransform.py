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

