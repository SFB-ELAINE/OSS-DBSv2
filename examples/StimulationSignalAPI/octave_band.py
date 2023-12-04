import matplotlib.pyplot as plt
import numpy as np

from ossdbs.stimulation_signals import (
    RectangleSignal,
    get_indices_in_octave_band,
    get_octave_band_indices,
    retrieve_time_domain_signal_from_fft,
)

frequency = 130
pulse_width = 600e-6
inter_pulse_width = 0.0
counter_pulse_width = 0.0

cutoff_frequency = 1e4

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)

# Full spectrum
timesteps, signal_retrieved = signal.retrieve_time_domain_signal(
    fft_signal, cutoff_frequency
)
plt.plot(timesteps / 1e-3, signal_retrieved, label="Full spectrum")

# Octave band

# only use positive frequencies
first_negative_freq = np.argwhere(fft_frequencies < 0)[0, 0]
frequencies = fft_frequencies[:first_negative_freq]
fourier_coefficients = fft_signal[:first_negative_freq]

frequency_indices = get_octave_band_indices(frequencies)
if not np.isclose(fourier_coefficients[0], 0.0):
    frequency_indices = np.insert(frequency_indices, 0, 0)

solution = np.zeros(fourier_coefficients.shape, dtype=complex)

for freq_idx in frequency_indices:
    band_indices = get_indices_in_octave_band(freq_idx, frequency_indices)
    print(band_indices)
    for oct_idx in band_indices:
        solution[oct_idx] = fourier_coefficients[oct_idx]

print(solution)

# Because we use full FFT we also need negative frequencies
n_frequencies = int(cutoff_frequency / frequency)
# to account for DC, too
tmp_freq_domain = np.zeros(n_frequencies + 1, dtype=complex)
if (n_frequencies + 1) % 2 == 1:  # if odd
    tmp_freq_domain = np.append(tmp_freq_domain, tmp_freq_domain[-1:0:-1])
else:
    tmp_freq_domain = np.append(tmp_freq_domain, tmp_freq_domain[-2:0:-1])

frequency_indices = frequencies / frequency
frequency_indices = frequency_indices.astype(np.uint16)
result_in_time = np.zeros(len(tmp_freq_domain))
# write frequencies
for idx, idx_freq in enumerate(frequency_indices):
    tmp_freq_domain[idx_freq] = solution[idx]
    # reverse order for negative frequencies
    if idx > 0:
        tmp_freq_domain[len(tmp_freq_domain) - idx_freq] = np.conjugate(solution[idx])
# convert to time domain
timesteps, tmp = retrieve_time_domain_signal_from_fft(
    tmp_freq_domain, cutoff_frequency, frequency
)
result_in_time[:] = tmp[:]

plt.plot(timesteps / 1e-3, result_in_time, label="Octave band")
plt.legend()
plt.show()
