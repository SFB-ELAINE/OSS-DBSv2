import json

import impedancefitter as ifit
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0 as e0

from ossdbs.dielectric_model import default_dielectric_parameters, dielectric_models

model = dielectric_models["ColeCole4"]
material_model = {}
for material, parameters in default_dielectric_parameters["ColeCole4"].items():
    material_model[material] = model(parameters)


def get_benchmark(tissue, frequencies):
    """Get benchmark from impedancefitter."""
    with open(tissue + "colecole4gabriel.json") as stream:
        resultcc4 = json.load(stream)
    omega = 2.0 * np.pi * frequencies
    ecmcc4 = ifit.get_equivalent_circuit_model("ColeCole4")
    Zcc4 = ecmcc4.eval(omega=omega, **resultcc4)
    c0 = resultcc4["c0"] * 1e-12
    eps_ast = 1.0 / (1j * 2.0 * np.pi * frequencies * c0 * Zcc4)
    eps_r, conductivity = ifit.utils.return_diel_properties(omega, Zcc4, c0)
    conductivity_ast = conductivity + 1j * omega * eps_r * e0
    return eps_ast * e0, conductivity_ast


frequencies = np.logspace(1, 6)
permittivity = np.zeros(len(frequencies), dtype=np.complex128)
conductivity = np.zeros(len(frequencies), dtype=np.complex128)
for idx, frequency in enumerate(frequencies):
    omega = 2 * np.pi * frequency
    permittivity[idx] = material_model["CSF"].complex_permittivity(omega)

for idx, frequency in enumerate(frequencies):
    omega = 2 * np.pi * frequency
    conductivity[idx] = material_model["CSF"].complex_conductivity(omega)

permittivity_bench, conductivity_bench = get_benchmark("csf", frequencies)
plt.xscale("log")
plt.plot(frequencies, conductivity.real, label="OSS-DBS")
plt.plot(frequencies, conductivity_bench.real, label="Benchmark", ls="dashed")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Real part of conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()

plt.xscale("log")
plt.plot(frequencies, conductivity.imag, label="OSS-DBS")
plt.plot(frequencies, conductivity_bench.imag, label="Benchmark", ls="dashed")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Imaginary part of conductivity / S m$^{-1}$")
plt.legend()
plt.tight_layout()
plt.show()

plt.xscale("log")
plt.yscale("log")
plt.plot(frequencies, permittivity.real / e0, label="OSS-DBS")
plt.plot(frequencies, permittivity_bench.real / e0, label="Benchmark", ls="dashed")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Relative permittivity")
plt.legend()
plt.tight_layout()
plt.show()

plt.xscale("log")
plt.plot(frequencies, permittivity.imag, label="OSS-DBS")
plt.plot(frequencies, permittivity_bench.imag, label="Benchmark", ls="dashed")
plt.xlabel("Frequency / Hz")
plt.ylabel(r"Imaginary part of permittivity$")
plt.legend()
plt.tight_layout()
plt.show()
