from ossdbs.stimulation_signals import RectangleSignal, generate_signal
import timeit

frequency = 130
pulse_width = 60e-6
inter_pulse_width = 120e-6
counter_pulse_width = 120e-6

cutoff_frequency = 1e5

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
frequencies, coefficients = signal.get_frequencies_and_fourier_coefficients(cutoff_frequency)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)

dt = 1.0 / cutoff_frequency
timesteps = int(cutoff_frequency / frequency)

print(timeit.timeit(lambda: generate_signal(coefficients, frequencies / frequency, frequency, dt, 0.0, timesteps, SigmaApprox=True), setup="pass", number=10))
print(timeit.timeit(lambda: signal.retrieve_time_domain_signal(fft_signal, cutoff_frequency), setup="pass", number=10))
