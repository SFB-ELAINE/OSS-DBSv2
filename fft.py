import numpy as np
import matplotlib.pyplot as plt
import time

n_samples = 10000
pulse_width = 0.0078
cut_off = 0.1
n_ifft = 1
time_step = 1 / n_samples

n_pulse_samples = int(n_samples * pulse_width)
print('n_pulse_samples', n_pulse_samples)
signal = np.zeros(n_samples)
signal[:n_pulse_samples] = 1

sp = np.fft.rfft(signal)
frequencies = np.fft.rfftfreq(len(signal), time_step / 130)
print(frequencies)
plt.plot(frequencies, sp.real, frequencies, sp.imag)
plt.show()

cut_off_samples = int(min(cut_off, 1) * n_samples / 2) + 1
s = np.fft.irfft(sp[:cut_off_samples]) * cut_off
print(len(s))
plt.plot(s)
plt.show()

sp2d = np.array([sp] * n_ifft)

t = time.time()
s_2d = np.fft.irfft(sp2d, axis=1)
duration = time.time() - t
print('duration', duration)
print(s_2d.shape)
