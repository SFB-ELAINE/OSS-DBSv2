from ossdbs.dielectric_model import DIELECTRIC_MODELS
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0 as e0

model = DIELECTRIC_MODELS["ColeCole4"]

frequencies = np.logspace(1, 6)
permittivity = np.zeros(len(frequencies), dtype=np.complex128)
conductivity = np.zeros(len(frequencies), dtype=np.complex128)
for idx, frequency in enumerate(frequencies):
    omega = 2 * np.pi * frequency
    permittivity[idx] = model.complex_permittivity("White matter", omega)

for idx, frequency in enumerate(frequencies):
    omega = 2 * np.pi * frequency
    conductivity[idx] = model.complex_conductivity("White matter", omega)

plt.xscale("log")
plt.plot(frequencies, conductivity.real)
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Real part of conductivity / S m$^{-1}$")
plt.show()

plt.xscale("log")
plt.plot(frequencies, conductivity.imag)
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Imaginary part of conductivity / S m$^{-1}$")
plt.show()

plt.xscale("log")
plt.yscale("log")
plt.plot(frequencies, permittivity.real / e0)
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Relative permittivity")
plt.show()

plt.xscale("log")
plt.plot(frequencies, permittivity.imag)
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Imaginary part of permittivity / S m$^{-1}$")
plt.show()
