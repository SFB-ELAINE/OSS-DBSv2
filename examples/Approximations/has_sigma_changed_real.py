import matplotlib.pyplot as plt
import numpy as np

from ossdbs.dielectric_model import default_dielectric_parameters, dielectric_models
from ossdbs.stimulation_signals import (
    RectangleSignal,
    get_octave_band_indices,
)
from ossdbs.utils import have_dielectric_properties_changed

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
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )
    conductivitiesCC3needed[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )
    conductivitiesCC3notneeded[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )

    conductivitiesCC4[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )
    conductivitiesCC4needed[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )
    conductivitiesCC4notneeded[material_test] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=float
    )

for idx, frequency in enumerate(frequencies):
    if idx > 0:
        changed_cc3 = have_dielectric_properties_changed(
            material_modelCC3,
            is_complex=False,
            old_freq=frequencies[idx - 1],
            new_freq=frequency,
            threshold=threshold,
        )
    else:
        changed_cc3 = True
    if idx > 0:
        changed_cc4 = have_dielectric_properties_changed(
            material_modelCC4,
            is_complex=False,
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
            ].conductivity(2.0 * np.pi * frequency)
        else:
            conductivitiesCC3notneeded[material_test][idx] = material_modelCC3[
                material_test
            ].conductivity(2.0 * np.pi * frequency)

        if changed_cc4:
            conductivitiesCC4needed[material_test][idx] = material_modelCC4[
                material_test
            ].conductivity(2.0 * np.pi * frequency)
        else:
            conductivitiesCC4notneeded[material_test][idx] = material_modelCC4[
                material_test
            ].conductivity(2.0 * np.pi * frequency)

        conductivitiesCC3[material_test][idx] = material_modelCC3[
            material_test
        ].conductivity(2.0 * np.pi * frequency)
        conductivitiesCC4[material_test][idx] = material_modelCC4[
            material_test
        ].conductivity(2.0 * np.pi * frequency)


for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(frequencies, conductivitiesCC3[material], label=material)
    plt.plot(frequencies, conductivitiesCC3needed[material], "o", color="green")
    plt.plot(frequencies, conductivitiesCC3notneeded[material], "o", color="red")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Real part of conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()

for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(frequencies, conductivitiesCC4[material], label=material)
    plt.plot(frequencies, conductivitiesCC4needed[material], "o", color="green")
    plt.plot(frequencies, conductivitiesCC4notneeded[material], "o", color="red")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Real part of conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()

plt.title("Compare real and imaginary part CC3")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(frequencies, conductivitiesCC3[material], label=material + " real part")
    plt.plot(frequencies, conductivitiesCC3needed[material], "o", color="green")
    plt.plot(frequencies, conductivitiesCC3notneeded[material], "o", color="red")

plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()


plt.title("Compare real and imaginary part CC4")
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(frequencies, conductivitiesCC4[material], label=material + " real part")
    plt.plot(frequencies, conductivitiesCC4needed[material], "o", color="green")
    plt.plot(frequencies, conductivitiesCC4notneeded[material], "o", color="red")

plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()

for material in material_modelCC4.keys():
    plt.plot(
        frequencies[1:],
        np.abs(np.diff(conductivitiesCC3[material]))
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
        np.abs(np.diff(conductivitiesCC4[material]))
        / np.abs(conductivitiesCC4[material][0:-1]),
        label=material,
    )
plt.legend()
plt.ylabel("Relative difference")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.show()
