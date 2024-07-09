import matplotlib.pyplot as plt
import numpy as np

from ossdbs.stimulation_signals import RectangleSignal, get_octave_band_indices

frequency = 130.0
pulse_width = 60e-6
cutoff_frequency = 1e6

signal = RectangleSignal(frequency, pulse_width, 0.0, 0.0)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)

# Select only the positive frequencies
# only use positive frequencies
first_negative_freq = np.argwhere(fft_frequencies < 0)[0, 0]
frequencies = fft_frequencies[:first_negative_freq]
print(fft_frequencies.shape)
fourier_coefficients = fft_signal[:first_negative_freq]
n_pos_frequencies = len(fourier_coefficients)
# even signal
signal_length = len(fft_signal)
if signal_length % 2 == 0:
    frequencies = np.append(
        frequencies, -1.0 * fft_frequencies[first_negative_freq + 1]
    )
    fourier_coefficients = np.append(
        fourier_coefficients, np.conjugate(fft_signal[first_negative_freq + 1])
    )

frequency_indices = get_octave_band_indices(frequencies)
print(frequency_indices)
print(frequencies[frequency_indices])

# Plot frequency domain signals
plt.stem(frequencies, np.abs(fourier_coefficients), markerfmt=" ")
plt.stem(
    frequencies[frequency_indices],
    np.abs(fourier_coefficients)[frequency_indices],
    markerfmt="C1o",
)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.xscale("log")
plt.title("FFT Spectrum")
plt.savefig("spectrum.png")
plt.show()
