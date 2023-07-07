from ossdbs.dielectric_models import DIELECTRIC_MODELS
import numpy as np
import matplotlib.pyplot as plt

model = DIELECTRIC_MODELS["ColeCole4"]

frequencies = np.logspace(1, 6)
permittivity = np.zeros(len(frequencies), dtype=np.complex)
conductivity = np.zeros(len(frequencies), dtype=np.complex)
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
