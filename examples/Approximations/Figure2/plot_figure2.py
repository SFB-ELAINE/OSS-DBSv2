"""Figure 2 — Monophasic stimulation pulse overview.

Four-panel figure (2x2):
  A: FFT coefficient amplitudes (stem plot, log x-axis)
  B: CC3 dielectric properties with frequency-skipping check (green=needed, red=skipped)
  C: CC3 ratio between real and imaginary conductivity with frequency-skipping check
  D: Time-domain signal reconstructed at different cutoff frequencies

Assumes a monophasic rectangular pulse, 60 us pulse width, 130 Hz.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0 as e0

from ossdbs.dielectric_model import default_dielectric_parameters, dielectric_models
from ossdbs.stimulation_signals import RectangleSignal, get_octave_band_indices
from ossdbs.stimulation_signals.utilities import adjust_cutoff_frequency
from ossdbs.utils import have_dielectric_properties_changed

# ---------- Signal parameters ----------
frequency = 130.0  # Hz
pulse_width = 60e-6  # s
inter_pulse_width = 0.0
counter_pulse_width = 0.0
cutoff_frequency = 1e6  # Hz
threshold = 0.01  # dielectric change threshold

signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)

# ---------- Figure setup ----------
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 10
plt.rcParams["axes.labelsize"] = 11
plt.rcParams["legend.fontsize"] = 8
plt.rcParams["text.latex.preamble"] = r"\usepackage{siunitx}"

fig, axes = plt.subplots(2, 2, figsize=(10, 7))
ax_a, ax_b, ax_c, ax_d = axes.flat

# ============================================================
# Panel A — FFT spectrum
# ============================================================
fft_frequencies, fft_signal, signal_length = signal.get_fft_spectrum(cutoff_frequency)
positive_indices = fft_frequencies >= 0
pos_freqs = fft_frequencies[positive_indices]
pos_coeffs = np.abs(fft_signal[positive_indices])

ax_a.stem(pos_freqs, pos_coeffs, markerfmt=" ", basefmt="C0-", linefmt="C0-")
ax_a.set_xlabel("Frequency / Hz")
ax_a.set_ylabel("FFT coefficient amplitude")
ax_a.set_xscale("log")
ax_a.set_title("A", loc="left", fontweight="bold")

# ============================================================
# Shared dielectric setup for panels B and C
# ============================================================
fft_frequencies_dc, _, _ = signal.get_fft_spectrum(cutoff_frequency)
frequency_indices = get_octave_band_indices(fft_frequencies_dc)
frequencies = frequency_indices * frequency

modelCC3 = dielectric_models["ColeCole3"]
modelCC4 = dielectric_models["ColeCole4"]
material_modelCC3 = {}
material_modelCC4 = {}
for material, parameters in default_dielectric_parameters["ColeCole3"].items():
    if material in ["Blood", "Unknown"]:
        continue
    material_modelCC3[material] = modelCC3(parameters)
    material_modelCC4[material] = modelCC4(
        default_dielectric_parameters["ColeCole4"][material]
    )


# Compute conductivities + classify needed/not-needed
def compute_dielectric_classification(material_models, frequencies, threshold):
    """Return full, needed, and not-needed complex conductivities per material."""
    full = {}
    needed = {}
    not_needed = {}
    for mat in material_models:
        full[mat] = np.full(frequencies.shape, np.nan, dtype=complex)
        needed[mat] = np.full(frequencies.shape, np.nan, dtype=complex)
        not_needed[mat] = np.full(frequencies.shape, np.nan, dtype=complex)

    for idx, freq in enumerate(frequencies):
        if idx > 0:
            changed = have_dielectric_properties_changed(
                material_models,
                is_complex=True,
                old_freq=frequencies[idx - 1],
                new_freq=freq,
                threshold=threshold,
            )
        else:
            changed = True

        for mat in material_models:
            sigma = material_models[mat].complex_conductivity(2.0 * np.pi * freq)
            full[mat][idx] = sigma
            if changed:
                needed[mat][idx] = sigma
            else:
                not_needed[mat][idx] = sigma

    return full, needed, not_needed


def plot_dielectric_panel(ax, material_models, frequencies, threshold, title):
    """Plot conductivity + permittivity with needed/skipped markers."""
    full, needed, not_needed = compute_dielectric_classification(
        material_models, frequencies, threshold
    )
    omega = 2.0 * np.pi * frequencies

    ax2 = ax.twinx()
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax.set_xscale("log")

    for mat in material_models:
        if mat == "CSF":
            continue
        (line,) = ax.plot(frequencies, full[mat].real, label=mat)
        ax.plot(frequencies, needed[mat].real, "o", color="green", markersize=3)
        ax.plot(frequencies, not_needed[mat].real, "o", color="red", markersize=3)
        perm = full[mat].imag / omega / e0
        perm_needed = needed[mat].imag / omega / e0
        perm_not_needed = not_needed[mat].imag / omega / e0
        ax2.plot(frequencies, perm, ls="dashed", color=line.get_color())
        ax2.plot(frequencies, perm_needed, "o", color="green", markersize=3)
        ax2.plot(frequencies, perm_not_needed, "o", color="red", markersize=3)

    ax.set_xlabel("Frequency / Hz")
    ax.set_ylabel(r"Real conductivity / S\,m$^{-1}$")
    ax2.set_ylabel("Rel. permittivity")
    ax.legend(loc="upper center", fontsize=7)
    ax.set_title(title, loc="left", fontweight="bold")


# ============================================================
# Panel B — CC3 dielectric properties
# ============================================================
plot_dielectric_panel(ax_b, material_modelCC3, frequencies, threshold, "B")

# ============================================================
# Panel C — Ratio between real and imaginary part (CC4)
# ============================================================
full_cc3, needed_cc3, not_needed_cc3 = compute_dielectric_classification(
    material_modelCC3, frequencies, threshold
)
ax_c.set_xscale("log")
ax_c.set_yscale("log")
for mat in material_modelCC4:
    if mat == "CSF":
        continue
    ratio = full_cc3[mat].real / np.abs(full_cc3[mat].imag)
    ratio_needed = needed_cc3[mat].real / np.abs(needed_cc3[mat].imag)
    ratio_not_needed = not_needed_cc3[mat].real / np.abs(not_needed_cc3[mat].imag)
    ax_c.plot(frequencies, ratio, label=mat)
    ax_c.plot(frequencies, ratio_needed, "o", color="green", markersize=3)
    ax_c.plot(frequencies, ratio_not_needed, "o", color="red", markersize=3)
ax_c.set_xlabel("Frequency / Hz")
ax_c.set_ylabel("Ratio between real and imaginary part")
ax_c.legend(loc="upper right", fontsize=7)
ax_c.set_title("C", loc="left", fontweight="bold")

# ============================================================
# Panel D — Time-domain signal at different cutoff frequencies
# ============================================================

# 1 MHz — solid line
adj_cutoff_1M = adjust_cutoff_frequency(2.0 * 1e6, frequency)
dt_1M = 1.0 / adj_cutoff_1M
n_steps_1M = int(adj_cutoff_1M / frequency)
td_1M = signal.get_time_domain_signal(dt_1M, n_steps_1M)
t_1M = dt_1M * np.arange(n_steps_1M)
ax_d.plot(
    t_1M * 1e6,
    td_1M,
    linestyle="-",
    linewidth=1.2,
    label=r"\SI{1}{\mega\hertz}",
)

# 100 kHz — dots only
adj_cutoff_100k = adjust_cutoff_frequency(2.0 * 1e5, frequency)
dt_100k = 1.0 / adj_cutoff_100k
n_steps_100k = int(adj_cutoff_100k / frequency)
td_100k = signal.get_time_domain_signal(dt_100k, n_steps_100k)
t_100k = dt_100k * np.arange(n_steps_100k)
ax_d.plot(
    t_100k * 1e6,
    td_100k,
    linestyle="none",
    marker="o",
    markersize=3,
    label=r"\SI{100}{\kilo\hertz}",
)
ax_d.set_xlabel(r"Time / $\mu$s")
ax_d.set_ylabel("Amplitude")
ax_d.set_xlim(left=0, right=5.0 * pulse_width * 1e6)
ax_d.legend(loc="upper right")
ax_d.set_title("D", loc="left", fontweight="bold")

# ---------- Save ----------
fig.tight_layout()
fig.savefig("figure2.pdf", bbox_inches="tight")
fig.savefig("figure2.svg", bbox_inches="tight")
plt.close()
print("Saved figure2.pdf / .svg")
