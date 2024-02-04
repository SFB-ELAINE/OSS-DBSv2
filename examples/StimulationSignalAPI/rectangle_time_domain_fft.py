from ossdbs.stimulation_signals import RectangleSignal
import matplotlib.pyplot as plt

frequency = 130
pulse_width = 60e-6
inter_pulse_width = 120e-6
counter_pulse_width = 120e-6

cutoff_frequency = 1e5

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)

timesteps, signal_retrieved = signal.retrieve_time_domain_signal(fft_signal, cutoff_frequency)
plt.plot(timesteps / 1e-3, signal_retrieved)
plt.show()
