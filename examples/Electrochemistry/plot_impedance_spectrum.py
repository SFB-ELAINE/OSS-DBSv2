"""Plot impedance spectrum of the Lempka 2009 electrode-tissue interface model.

The model from Lempka et al. (2009) "In vivo impedance spectroscopy of deep
brain stimulation electrodes" is:

    Z = CPE_dl + R_inf + parallel(R_delta, CPE_tissue)

where CPE_dl is the double-layer constant-phase element, R_inf is the
high-frequency resistance, and parallel(R_delta, CPE_tissue) captures the
tissue relaxation.

Parameter sets correspond to different time points post-implantation
(day 1, day 15) and pre/post stimulation conditions.

Usage
-----
    python plot_impedance_spectrum.py
"""

import impedancefitter as ifit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Frequency range (Hz) ---
f_min = 10.0
f_max = 100e3
n_points = 200


# --- Lempka 2009 equivalent circuit models ---
lempka_cpe = "CPE_dl"
lempka_tissue = "R_inf + parallel(R_delta, CPE_tissue)"
lempka_full = lempka_cpe + "+ " + lempka_tissue

ecm_full = ifit.get_equivalent_circuit_model(lempka_full)
ecm_tissue = ifit.get_equivalent_circuit_model(lempka_tissue)
ecm_cpe = ifit.get_equivalent_circuit_model(lempka_cpe)

# Load results

sim_impedance = pd.read_csv("Results_LempkaImpedance_spectrum/impedance.csv")
sim_frequencies = sim_impedance["freq"].to_numpy()
omega = 2.0 * np.pi * sim_frequencies
sim_Z = sim_impedance["real"].to_numpy() + 1j * sim_impedance["imag"].to_numpy()
sim_impedance_no_interface = pd.read_csv(
    "Results_NoInterface_Lempka_Spectrum/impedance.csv"
)
assert np.all(
    np.isclose(sim_frequencies, sim_impedance_no_interface["freq"].to_numpy())
)
sim_no_interface_Z = (
    sim_impedance_no_interface["real"].to_numpy()
    + 1j * sim_impedance_no_interface["imag"].to_numpy()
)

# --- Parameter sets from Lempka 2009 ---
# Note: R_delta and CPE_tissue parameters are rough estimates from the paper
param_sets = {
    "Day 1": {
        "dl_k": 1.5e6,
        "dl_alpha": 0.8,
        "inf_R": 1.25e3,
        "delta_R": 20e3 / 11.0,
        "tissue_k": 10e6 / 67.0,
        "tissue_alpha": 0.5,
    },
    "Day 15": {
        "dl_k": 1e6,
        "dl_alpha": 0.75,
        "inf_R": 1.35e3,
        "delta_R": 20e3,
        "tissue_k": 10e6,
        "tissue_alpha": 0.75,
    },
    "Before stim": {
        "dl_k": 1e6,
        "dl_alpha": 0.75,
        "inf_R": 2.5e3,
        "delta_R": 28e3,
        "tissue_k": 18e6,
        "tissue_alpha": 0.75,
    },
    "After stim": {
        "dl_k": 1.2e6,
        "dl_alpha": 0.8,
        "inf_R": 2.0e3,
        "delta_R": 5e3,
        "tissue_k": 10e6,
        "tissue_alpha": 0.5,
    },
}

# --- Compute impedance for each parameter set ---
Z_full = {}
Z_tissue = {}
Z_cpe = {}
for name, params in param_sets.items():
    Z_full[name] = ecm_full.eval(omega=omega, **params)
    Z_tissue[name] = ecm_tissue.eval(omega=omega, **params)
    Z_cpe[name] = ecm_cpe.eval(omega=omega, **params)

# --- Bode plot: |Z| and phase ---
fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

for name in param_sets:
    Z = Z_full[name]
    ax_mag.plot(sim_frequencies, np.abs(Z), label=name)
    ax_phase.plot(sim_frequencies, np.angle(Z, deg=True), label=name)
ax_mag.plot(sim_frequencies, np.abs(sim_Z), label="Simulation")
ax_phase.plot(sim_frequencies, np.angle(sim_Z, deg=True), label="Simulation")

ax_mag.set_ylabel(r"|Z| / $\Omega$")
ax_mag.set_xscale("log")
ax_mag.set_yscale("log")
ax_mag.legend()
ax_mag.grid(True, which="both", ls=":", alpha=0.5)
ax_mag.set_title("Lempka 2009 — Electrode-tissue interface impedance")

ax_phase.set_xlabel("Frequency / Hz")
ax_phase.set_ylabel("Phase / deg")
ax_phase.set_xscale("log")
ax_phase.grid(True, which="both", ls=":", alpha=0.5)

fig.tight_layout()
fig.savefig("impedance_spectrum_bode.pdf")
fig.savefig("impedance_spectrum_bode.svg")
plt.show()

# --- Tissue vs full model comparison (Day 1) ---
fig2, (ax2_mag, ax2_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

ref_name = "Day 1"
ax2_mag.plot(sim_frequencies, np.abs(sim_Z), label="Full model")
ax2_mag.plot(
    sim_frequencies,
    np.abs(sim_no_interface_Z + Z_cpe[ref_name]),
    label="Full model (estimated)",
)
ax2_mag.plot(
    sim_frequencies, np.abs(sim_no_interface_Z), label="Tissue only", ls="dashed"
)
ax2_mag.set_ylabel(r"|Z| / $\Omega$")
ax2_mag.set_xscale("log")
ax2_mag.set_yscale("log")
ax2_mag.legend()
ax2_mag.grid(True, which="both", ls=":", alpha=0.5)
ax2_mag.set_title(f"Lempka 2009 — Full model vs tissue contribution ({ref_name})")

ax2_phase.plot(sim_frequencies, np.angle(sim_Z, deg=True), label="Full model")
ax2_phase.plot(
    sim_frequencies,
    np.angle(sim_no_interface_Z + Z_cpe[ref_name], deg=True),
    label="Full model (estimated)",
)
ax2_phase.plot(
    sim_frequencies,
    np.angle(sim_no_interface_Z, deg=True),
    label="Tissue only",
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
    # Z = Z_full[name]
    Z = sim_no_interface_Z + Z_cpe[name]
    ax3.plot(Z.real, -Z.imag, label=name)
ax3.plot(sim_Z.real, -sim_Z.imag, label="Simulation")

ax3.set_xlabel(r"Re(Z) / $\Omega$")
ax3.set_ylabel(r"$-$Im(Z) / $\Omega$")
ax3.set_aspect("equal")
ax3.legend()
ax3.grid(True, ls=":", alpha=0.5)
ax3.set_title("Lempka 2009 — Nyquist plot")

fig3.tight_layout()
fig3.savefig("impedance_spectrum_nyquist.pdf")
fig3.savefig("impedance_spectrum_nyquist.svg")
plt.show()
