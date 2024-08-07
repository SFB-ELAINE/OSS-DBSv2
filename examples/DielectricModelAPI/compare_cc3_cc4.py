import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0 as e0

from ossdbs.dielectric_model import default_dielectric_parameters, dielectric_models

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
conductivitiesCC4 = {}

frequencies = np.logspace(1, 7)

for material in material_modelCC4.keys():
    conductivitiesCC3[material] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )
    conductivitiesCC4[material] = np.full(
        shape=frequencies.shape, fill_value=np.nan, dtype=complex
    )

    for idx, frequency in enumerate(frequencies):
        conductivitiesCC3[material][idx] = material_modelCC3[
            material
        ].complex_conductivity(2.0 * np.pi * frequency)
        conductivitiesCC4[material][idx] = material_modelCC4[
            material
        ].complex_conductivity(2.0 * np.pi * frequency)

omega = 2.0 * np.pi * frequencies
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
counter = 0
ax1.set_yscale("log")
ax2.set_yscale("log")
plt.xscale("log")

for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    ax1.plot(frequencies, conductivitiesCC3[material].real, c=colors[counter * 2])

    ax1.plot(
        frequencies,
        conductivitiesCC4[material].real,
        "v-",
        markevery=5,
        c=colors[counter * 2],
    )
    # plot on second axis
    ax2.plot(
        frequencies,
        conductivitiesCC3[material].imag / omega / e0,
        "-.",
        color=colors[counter * 2],
    )
    ax2.plot(
        frequencies,
        conductivitiesCC4[material].imag / omega / e0,
        ls="-.",
        marker="o",
        markevery=5,
        c=colors[counter * 2],
        label=material,
    )
    counter += 1
ax1.set_xlabel("Frequency / Hz")
ax1.set_ylabel(r"Conductivity / S m$^{-1}$")
ax2.set_ylabel("Rel. permittivity")
plt.legend(loc="upper center")
plt.tight_layout()
plt.savefig("CC3_vs_CC4.pdf")
# plt.show()
plt.close()

counter = 0
for material in material_modelCC4.keys():
    if material == "CSF":
        continue
    # Plotting conductivity
    plt.plot(
        frequencies,
        conductivitiesCC3[material].real / np.abs(conductivitiesCC3[material].imag),
        c=colors[counter * 2],
    )
    plt.plot(
        frequencies,
        conductivitiesCC4[material].real / np.abs(conductivitiesCC4[material].imag),
        "o-",
        markevery=5,
        label=material,
        c=colors[counter * 2],
    )
    counter += 1
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency / Hz")
plt.ylabel("Ratio between real and imaginary part")
plt.legend()
plt.tight_layout()
plt.savefig("ratio_CC3_vs_CC4.pdf")
plt.show()
