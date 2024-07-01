import matplotlib.pyplot as plt
import numpy as np

from ossdbs.stimulation_signals import RectangleSignal, TrapezoidSignal, TriangleSignal

settings = {
    "StimulationSignal": {
        "Type": "Rectangle",
        "Frequency[Hz]": 130.0,
        "PulseWidth[us]": 60e-6,
        "CounterPulseWidth[us]": 120e-6,
        "InterPulseWidth[us]": 120e-6,
        "CutoffFrequency": 1e5,
        "TopWidth[us]": 30e-6,
    }
}

signal_type = settings["StimulationSignal"]["Type"]
frequency = settings["StimulationSignal"]["Frequency[Hz]"]
pulse_width = settings["StimulationSignal"]["PulseWidth[us]"]
counter_pulse_width = settings["StimulationSignal"]["CounterPulseWidth[us]"]
inter_pulse_width = settings["StimulationSignal"]["InterPulseWidth[us]"]
cutoff_frequency = settings["StimulationSignal"]["CutoffFrequency"]
top_width = settings["StimulationSignal"]["TopWidth[us]"]

dt = 1.0 / cutoff_frequency
timesteps = int(cutoff_frequency / frequency)

if signal_type == "Rectangle":
    signal = RectangleSignal(
        frequency, pulse_width, inter_pulse_width, counter_pulse_width
    )
elif signal_type == "Triangle":
    signal = TriangleSignal(
        frequency, pulse_width, inter_pulse_width, counter_pulse_width
    )
elif signal_type == "Trapezoid":
    signal = TrapezoidSignal(
        frequency, pulse_width, inter_pulse_width, top_width, counter_pulse_width
    )

time_domain_signal = signal.get_time_domain_signal(dt, timesteps)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)

# Select only the positive frequencies
positive_indices = fft_frequencies >= 0
positive_frequencies = fft_frequencies[positive_indices]
positive_fft_signal = np.abs(fft_signal[positive_indices])

# Plot frequency domain signals
plt.stem(positive_frequencies, positive_fft_signal, markerfmt=" ")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.xscale("log")
plt.title("FFT Spectrum")
plt.show()

# Plot time domain signals (original and retrieved ones)
timesteps_retrived, retrived_signal = signal.retrieve_time_domain_signal(
    fft_signal, cutoff_frequency
)
plt.plot(timesteps_retrived, retrived_signal, label="retrieved signal")
plt.plot(dt * np.arange(timesteps), time_domain_signal, label="original signal")
plt.legend()
plt.show()
