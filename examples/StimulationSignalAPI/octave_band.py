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
cutoff_frequency = 10010

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)
print(fft_frequencies.shape)
print(fft_signal.shape)
signal_length = len(fft_signal)

# Full spectrum
timesteps, signal_retrieved = signal.retrieve_time_domain_signal(
    fft_signal, cutoff_frequency
)
plt.plot(timesteps / 1e-3, signal_retrieved, label="Full spectrum ideal")


# only use positive frequencies
first_negative_freq = np.argwhere(fft_frequencies < 0)[0, 0]
frequencies = fft_frequencies[:first_negative_freq]
print(fft_frequencies.shape)
fourier_coefficients = fft_signal[:first_negative_freq]
n_pos_frequencies = len(fourier_coefficients)
# even signal
if signal_length % 2 == 0:
    frequencies = np.append(
        frequencies, -1.0 * fft_frequencies[first_negative_freq + 1]
    )
    fourier_coefficients = np.append(
        fourier_coefficients, np.conjugate(fft_signal[first_negative_freq + 1])
    )

frequency_indices = frequencies / frequency
frequency_indices = frequency_indices.astype(np.uint16)

tmp_freq_domain = np.zeros(signal_length, dtype=complex)

# write frequencies
for idx, idx_freq in enumerate(frequency_indices):
    # write frequencies but skip last frequency for odd signal
    if idx_freq == n_pos_frequencies:
        if signal_length % 2 == 1:
            tmp_freq_domain[idx_freq] = fourier_coefficients[idx]
    else:
        tmp_freq_domain[idx_freq] = fourier_coefficients[idx]
    # reverse order for negative frequencies
    if idx > 0:
        tmp_freq_domain[len(tmp_freq_domain) - idx_freq] = np.conjugate(
            fourier_coefficients[idx]
        )
# convert to time domain
timesteps, result_in_time = retrieve_time_domain_signal_from_fft(
    tmp_freq_domain, cutoff_frequency, frequency
)

plt.plot(timesteps / 1e-3, result_in_time, ls="-", label="Full spectrum")

# Octave band

frequency_indices = get_octave_band_indices(frequencies)
if not np.isclose(fourier_coefficients[0], 0.0):
    frequency_indices = np.insert(frequency_indices, 0, 0)

print(frequency_indices)

solution = np.zeros(signal_length, dtype=complex)

for freq_idx in frequency_indices:
    band_indices = get_indices_in_octave_band(
        freq_idx, frequency_indices, int(cutoff_frequency / frequency)
    )
    for oct_idx in band_indices:
        solution[oct_idx] = fourier_coefficients[oct_idx]
    print(f"Highest frequency in band: {frequencies[max(band_indices)]}")

# to account for DC, too
tmp_freq_domain = np.zeros(signal_length, dtype=complex)

frequency_indices = frequencies / frequency
frequency_indices = frequency_indices.astype(np.uint16)
print(frequency_indices)
result_in_time = np.zeros(len(tmp_freq_domain))

# write frequencies
for idx, idx_freq in enumerate(frequency_indices):
    # write frequencies but skip last frequency for odd signal
    if idx_freq == n_pos_frequencies:
        if signal_length % 2 == 1:
            tmp_freq_domain[idx_freq] = fourier_coefficients[idx]
    else:
        tmp_freq_domain[idx_freq] = fourier_coefficients[idx]
    # reverse order for negative frequencies
    if idx > 0:
        tmp_freq_domain[len(tmp_freq_domain) - idx_freq] = np.conjugate(
            fourier_coefficients[idx]
        )
# convert to time domain
timesteps, result_in_time = retrieve_time_domain_signal_from_fft(
    tmp_freq_domain, cutoff_frequency, frequency
)

plt.plot(timesteps / 1e-3, result_in_time, ls="dotted", label="Octave band")
plt.legend()
plt.show()
