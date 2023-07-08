from ossdbs.stimulation_signals import RectangleSignal
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
time_domain_signal = signal.get_time_domain_signal(dt, timesteps)
plt.plot(dt * np.arange(timesteps) / 1e-6, time_domain_signal)
plt.show()
