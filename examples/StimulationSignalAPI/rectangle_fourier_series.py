from ossdbs.stimulation_signals import RectangleSignal, generate_signal
import matplotlib.pyplot as plt
import numpy as np

frequency = 130
pulse_width = 60e-6
inter_pulse_width = 120e-6
counter_pulse_width = 120e-6

cutoff_frequency = 1e5
dt = 2e-6
timesteps = int(1. / frequency / dt)

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
frequencies, coefficients = signal.get_frequencies_and_fourier_coefficients(cutoff_frequency)

signal = generate_signal(coefficients, frequencies / frequency, frequency, dt, 0.0, timesteps, SigmaApprox=True)
signal_ringing = generate_signal(coefficients, frequencies / frequency, frequency, dt, 0.0, timesteps, SigmaApprox=False)
plt.xlim(-10, 500)
plt.plot(dt * np.arange(len(signal)) / 1e-6, signal_ringing)
plt.plot(dt * np.arange(len(signal)) / 1e-6, signal)
plt.show()
