import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0 as e0

from ossdbs.dielectric_model import default_dielectric_parameters, dielectric_models
from ossdbs.stimulation_signals import (
    RectangleSignal,
    get_octave_band_indices,
)
from ossdbs.utils import have_dielectric_properties_changed

# whether to show plots
show = False
frequency = 130
pulse_width = 60e-6
inter_pulse_width = 0.0
counter_pulse_width = 0.0

cutoff_frequency = 1e6
threshold = 0.01

# get frequencies
signal = RectangleSignal(frequency, pulse_width, inter_pulse_width, counter_pulse_width)
fft_frequencies, fft_signal = signal.get_fft_spectrum(cutoff_frequency)
first_negative_freq = np.argwhere(fft_frequencies < 0)[0, 0]
frequencies = fft_frequencies[:first_negative_freq]
frequency_indices = get_octave_band_indices(frequencies)
frequencies = frequency_indices * frequency

modelCC4 = dielectric_models["ColeCole4"]
modelCC3 = dielectric_models["ColeCole3"]
material_modelCC4 = {}
material_modelCC3 = {}
for material, parameters in default_dielectric_parameters["ColeCole3"].items():
    # we exclude blood
    if material in ["Blood", "Unknown"]:
        continue
    material_modelCC3[material] = modelCC3(parameters)
    material_modelCC4[material] = modelCC4(
        default_dielectric_parameters["ColeCole4"][material]
    )


conductivitiesCC3 = {}
conductivitiesCC3needed = {}
conductivitiesCC3notneeded = {}

conductivitiesCC4 = {}
conductivitiesCC4needed = {}
conductivitiesCC4notneeded = {}


for material_test in material_modelCC4.keys():
    conductivitiesCC3[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )
    conductivitiesCC3needed[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )
    conductivitiesCC3notneeded[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )

    conductivitiesCC4[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )
    conductivitiesCC4needed[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )
    conductivitiesCC4notneeded[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )

for idx, frequency in enumerate(frequencies):
    if idx > 0:
        changed_cc3 = have_dielectric_properties_changed(
            material_modelCC3,
            is_complex=True,
            old_freq=frequencies[idx - 1],
            new_freq=frequency,
            threshold=threshold,
        )
    else:
        changed_cc3 = True
    if idx > 0:
        changed_cc4 = have_dielectric_properties_changed(
            material_modelCC4,
            is_complex=True,
            old_freq=frequencies[idx - 1],
            new_freq=frequency,
            threshold=threshold,
        )
    else:
        changed_cc4 = True

    for material_test in material_modelCC4.keys():
        if changed_cc3:
            conductivitiesCC3needed[material_test][idx] = material_modelCC3[
                material_test
            ].complex_conductivity(2.0 * np.pi * frequency)
        else:
            conductivitiesCC3notneeded[material_test][idx] = material_modelCC3[
                material_test
            ].complex_conductivity(2.0 * np.pi * frequency)

        if changed_cc4:
            conductivitiesCC4needed[material_test][idx] = material_modelCC4[
                material_test
            ].complex_conductivity(2.0 * np.pi * frequency)
        else:
            conductivitiesCC4notneeded[material_test][idx] = material_modelCC4[
                material_test
            ].complex_conductivity(2.0 * np.pi * frequency)

        conductivitiesCC3[material_test][idx] = material_modelCC3[
            material_test
        ].complex_conductivity(2.0 * np.pi * frequency)
        conductivitiesCC4[material_test][idx] = material_modelCC4[
            material_test
        ].complex_conductivity(2.0 * np.pi * frequency)


omega = 2.0 * np.pi * frequencies
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.xscale("log")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    ax1.plot(frequencies, conductivitiesCC3[material].real, label=material)
    ax1.plot(frequencies, conductivitiesCC3needed[material].real, "o", color="green")
    ax1.plot(frequencies, conductivitiesCC3notneeded[material].real, "o", color="red")
    ax2.plot(frequencies, conductivitiesCC3[material].imag / omega / e0, ls="dashed")
    ax2.plot(
        frequencies,
        conductivitiesCC3needed[material].imag / omega / e0,
        "o",
        color="green",
    )
    ax2.plot(
        frequencies,
        conductivitiesCC3notneeded[material].imag / omega / e0,
        "o",
        color="red",
    )

plt.xscale("log")
ax1.set_ylabel(r"Real part of conductivity / S m$^{-1}$")
ax2.set_ylabel("Rel. permittivity")
ax1.legend(loc="upper center")
ax1.set_xlabel("Frequency / Hz")
plt.tight_layout()
plt.savefig("CC3_dielectric_properties.pdf")
if show:
    plt.show()
else:
    plt.close()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.xscale("log")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    ax1.plot(frequencies, conductivitiesCC4[material].real, label=material)
    ax1.plot(frequencies, conductivitiesCC4needed[material].real, "o", color="green")
    ax1.plot(frequencies, conductivitiesCC4notneeded[material].real, "o", color="red")
    # Plot permittivity
    ax2.plot(frequencies, conductivitiesCC4[material].imag / omega / e0, ls="dashed")
    ax2.plot(
        frequencies,
        conductivitiesCC4needed[material].imag / omega / e0,
        "o",
        color="green",
    )
    ax2.plot(
        frequencies,
        conductivitiesCC4notneeded[material].imag / omega / e0,
        "o",
        color="red",
    )

plt.yscale("log")
plt.xscale("log")
ax1.set_ylabel(r"Real part of conductivity / S m$^{-1}$")
ax2.set_ylabel("Rel. permittivity")
ax1.legend(loc="upper center")
ax1.set_xlabel("Frequency / Hz")
plt.tight_layout()
plt.savefig("CC4_dielectric_properties.pdf")
if show:
    plt.show()
else:
    plt.close()

print("Compare real and imaginary part CC3")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(
        frequencies,
        conductivitiesCC3[material].real / np.abs(conductivitiesCC3[material].imag),
        label=material,
    )
    plt.plot(
        frequencies,
        conductivitiesCC3needed[material].real
        / np.abs(conductivitiesCC3needed[material].imag),
        "o",
        color="green",
    )
    plt.plot(
        frequencies,
        conductivitiesCC3notneeded[material].real
        / np.abs(conductivitiesCC3notneeded[material].imag),
        "o",
        color="red",
    )
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel("Ratio between real and imaginary part")
plt.legend()
plt.tight_layout()
plt.savefig("CC3_ratios_real_imag.pdf")
if show:
    plt.show()
else:
    plt.close()


print("Compare real and imaginary part CC4")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(
        frequencies,
        conductivitiesCC4[material].real / np.abs(conductivitiesCC4[material].imag),
        label=material + " real part",
    )
    plt.plot(
        frequencies,
        conductivitiesCC4needed[material].real
        / np.abs(conductivitiesCC4needed[material].imag),
        "o",
        color="green",
    )
    plt.plot(
        frequencies,
        conductivitiesCC4notneeded[material].real
        / np.abs(conductivitiesCC4notneeded[material].imag),
        "o",
        color="red",
    )
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel("Ratio between real and imaginary part")
plt.legend()
plt.tight_layout()
plt.savefig("CC4_ratios_real_imag.pdf")
if show:
    plt.show()
else:
    plt.close()

"""
for material in material_modelCC4.keys():
    plt.plot(
        frequencies[1:],
        np.abs(np.diff(conductivitiesCC3[material].real))
        / np.abs(conductivitiesCC3[material][0:-1]),
        label=material,
    )
plt.legend()
plt.ylabel("Relative difference")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.show()

for material in material_modelCC4.keys():
    plt.plot(
        frequencies[1:],
        np.abs(np.diff(conductivitiesCC4[material].real))
        / np.abs(conductivitiesCC4[material][0:-1]),
        label=material,
    )
plt.legend()
plt.ylabel("Relative difference")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.show()
"""
