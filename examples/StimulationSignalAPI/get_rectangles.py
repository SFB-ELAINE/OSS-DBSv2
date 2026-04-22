import matplotlib.pyplot as plt
import numpy as np

from ossdbs.stimulation_signals import RectangleSignal
from ossdbs.stimulation_signals.utilities import adjust_cutoff_frequency

settings = {
    "StimulationSignal": {
        "Type": "Rectangle",
        "Frequency[Hz]": 130.0,
        "PulseWidth[us]": 60e-6,
        "CounterPulseWidth[us]": 0e-6,
        "InterPulseWidth[us]": 0e-6,
        "CutoffFrequency": 1e6,
        "TopWidth[us]": 0e-6,
    }
}

plt.rcParams["font.size"] = 16
plt.rcParams["legend.fontsize"] = 14
plt.rcParams["legend.framealpha"] = 1.0

signal_type = settings["StimulationSignal"]["Type"]
frequency = settings["StimulationSignal"]["Frequency[Hz]"]
pulse_width = settings["StimulationSignal"]["PulseWidth[us]"]
counter_pulse_width = settings["StimulationSignal"]["CounterPulseWidth[us]"]
inter_pulse_width = settings["StimulationSignal"]["InterPulseWidth[us]"]
cutoff_frequency = settings["StimulationSignal"]["CutoffFrequency"]
top_width = settings["StimulationSignal"]["TopWidth[us]"]

adj_cutoff_frequency = adjust_cutoff_frequency(2.0 * cutoff_frequency, frequency)
dt = 1.0 / adj_cutoff_frequency
timesteps = int(adj_cutoff_frequency / frequency)

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)

time_domain_signal = signal.get_time_domain_signal(dt, timesteps)
plt.plot(1e6 * dt * np.arange(timesteps), time_domain_signal, label="Original signal")
plt.xlim(left=0, right=1e6 * (10 * pulse_width + counter_pulse_width))
plt.xlabel(r"Time / $\mu$s")
plt.ylabel("Amplitude / V")
plt.tight_layout()
plt.savefig("time_domain_signal_mono.pdf")
plt.show()

# make biphasic
counter_pulse_width = 60e-6
signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
time_domain_signal = signal.get_time_domain_signal(dt, timesteps)
plt.plot(1e6 * dt * np.arange(timesteps), time_domain_signal, label="Original signal")
plt.xlim(left=0, right=1e6 * (10 * pulse_width + counter_pulse_width))
plt.xlabel(r"Time / $\mu$s")
plt.ylabel("Amplitude / V")
plt.tight_layout()
plt.savefig("time_domain_signal_bi.pdf")
plt.show()
