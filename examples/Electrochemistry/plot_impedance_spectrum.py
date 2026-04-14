"""Compare FEM simulation against tissue-FEM + CPE decomposition.

The reference setup follows Lempka et al. (2009), "In vivo impedance
spectroscopy of deep brain stimulation electrodes", in which the
macroscopic impedance measured between two DBS contacts is decomposed
into an electrode-tissue double-layer term and a tissue term:

    Z_measured ≈ Z_CPE_dl + Z_tissue

Here the tissue impedance is taken directly from a FEM simulation
without surface impedance (pure Dirichlet BC), so the decomposition
becomes

    Z_model(f) = Z_tissue_FEM(f) + Z_CPE_dl(f; dl_k, dl_alpha)

which is then compared against the full FEM simulation with the
``CPE_dl`` Robin BC on the active contact. Different ``(dl_k, dl_alpha)``
parameter sets from Lempka 2009 are plotted to illustrate the sensitivity
of the decomposition to the double-layer parameters.

Inputs (both computed over the same frequency range):

- ``Results_LempkaImpedance_spectrum/impedance.csv`` — full FEM with
  ``CPE_dl`` (``dl_k=1.5e6, dl_alpha=0.8``, i.e. the Day 1 set).
- ``Results_NoInterface_Lempka_Spectrum/impedance.csv`` — pure Dirichlet,
  same geometry/tissue/frequencies.

Usage
-----
    python plot_impedance_spectrum.py
"""

import impedancefitter as ifit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Lempka 2009 double-layer CPE only (the tissue part is taken from FEM) ---
ecm_cpe = ifit.get_equivalent_circuit_model("CPE_dl")

# --- Load FEM results ---
sim_impedance = pd.read_csv("Results_LempkaImpedance_spectrum/impedance.csv")
sim_frequencies = sim_impedance["freq"].to_numpy()
omega = 2.0 * np.pi * sim_frequencies
sim_Z = sim_impedance["real"].to_numpy() + 1j * sim_impedance["imag"].to_numpy()

sim_impedance_no_interface = pd.read_csv(
    "Results_NoInterface_Lempka_Spectrum/impedance.csv"
)
assert np.all(
    np.isclose(sim_frequencies, sim_impedance_no_interface["freq"].to_numpy())
), "Frequency grids of the two simulations must match."
sim_no_interface_Z = (
    sim_impedance_no_interface["real"].to_numpy()
    + 1j * sim_impedance_no_interface["imag"].to_numpy()
)

# --- Lempka 2009 CPE_dl parameter sets (dl_k, dl_alpha only) ---
# Only Day 1 matches the parameters used in the full FEM simulation.
param_sets = {
    "Day 1": {"dl_k": 1.5e6, "dl_alpha": 0.8},
    "Day 15": {"dl_k": 1e6, "dl_alpha": 0.75},
    "Before stim": {"dl_k": 1e6, "dl_alpha": 0.75},
    "After stim": {"dl_k": 1.2e6, "dl_alpha": 0.8},
}

# Reconstructed impedance: FEM tissue + ECM CPE_dl.
Z_cpe = {name: ecm_cpe.eval(omega=omega, **p) for name, p in param_sets.items()}
Z_model = {name: sim_no_interface_Z + Z_cpe[name] for name in param_sets}

# --- Bode plot: full FEM vs FEM-tissue + CPE_dl (all parameter sets) ---
fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

for name in param_sets:
    ax_mag.plot(sim_frequencies, np.abs(Z_model[name]), label=name)
    ax_phase.plot(sim_frequencies, np.angle(Z_model[name], deg=True), label=name)
ax_mag.plot(sim_frequencies, np.abs(sim_Z), label="Full FEM", color="k", ls="dashed")
ax_phase.plot(
    sim_frequencies, np.angle(sim_Z, deg=True), label="Full FEM", color="k", ls="dashed"
)

ax_mag.set_ylabel(r"|Z| / $\Omega$")
ax_mag.set_xscale("log")
ax_mag.set_yscale("log")
ax_mag.legend()
ax_mag.grid(True, which="both", ls=":", alpha=0.5)
ax_mag.set_title("Full FEM vs FEM-tissue + CPE_dl (Lempka 2009 parameter sets)")

ax_phase.set_xlabel("Frequency / Hz")
ax_phase.set_ylabel("Phase / deg")
ax_phase.set_xscale("log")
ax_phase.grid(True, which="both", ls=":", alpha=0.5)

fig.tight_layout()
fig.savefig("impedance_spectrum_bode.pdf")
fig.savefig("impedance_spectrum_bode.svg")
plt.show()

# --- Decomposition plot for the matching parameter set (Day 1) ---
fig2, (ax2_mag, ax2_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

ref_name = "Day 1"
ax2_mag.plot(sim_frequencies, np.abs(sim_Z), label="Full FEM")
ax2_mag.plot(
    sim_frequencies,
    np.abs(Z_model[ref_name]),
    label=f"FEM tissue + CPE_dl ({ref_name})",
)
ax2_mag.plot(
    sim_frequencies, np.abs(sim_no_interface_Z), label="FEM tissue only", ls="dashed"
)
ax2_mag.set_ylabel(r"|Z| / $\Omega$")
ax2_mag.set_xscale("log")
ax2_mag.set_yscale("log")
ax2_mag.legend()
ax2_mag.grid(True, which="both", ls=":", alpha=0.5)
ax2_mag.set_title(f"Impedance decomposition ({ref_name} CPE_dl parameters)")

ax2_phase.plot(sim_frequencies, np.angle(sim_Z, deg=True), label="Full FEM")
ax2_phase.plot(
    sim_frequencies,
    np.angle(Z_model[ref_name], deg=True),
    label=f"FEM tissue + CPE_dl ({ref_name})",
)
ax2_phase.plot(
    sim_frequencies,
    np.angle(sim_no_interface_Z, deg=True),
    label="FEM tissue only",
    ls="dashed",
)
ax2_phase.set_xlabel("Frequency / Hz")
ax2_phase.set_ylabel("Phase / deg")
ax2_phase.set_xscale("log")
ax2_phase.legend()
ax2_phase.grid(True, which="both", ls=":", alpha=0.5)

fig2.tight_layout()
fig2.savefig("impedance_spectrum_decomposition.pdf")
fig2.savefig("impedance_spectrum_decomposition.svg")
plt.show()

# --- Nyquist plot ---
fig3, ax3 = plt.subplots(figsize=(7, 6))

for name in param_sets:
    ax3.plot(Z_model[name].real, -Z_model[name].imag, label=name)
ax3.plot(sim_Z.real, -sim_Z.imag, label="Full FEM", color="k", ls="dashed")

ax3.set_xlabel(r"Re(Z) / $\Omega$")
ax3.set_ylabel(r"$-$Im(Z) / $\Omega$")
ax3.set_aspect("equal")
ax3.legend()
ax3.grid(True, ls=":", alpha=0.5)
ax3.set_title("Nyquist: full FEM vs FEM-tissue + CPE_dl")

fig3.tight_layout()
fig3.savefig("impedance_spectrum_nyquist.pdf")
fig3.savefig("impedance_spectrum_nyquist.svg")
plt.show()
