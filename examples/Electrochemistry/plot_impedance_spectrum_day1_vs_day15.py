"""Compare Day 1 vs Day 15 Lempka 2009 CPE_dl spectra.

Two full FEM simulations with different double-layer CPE parameters
from Lempka 2009 are compared against the decomposition

    Z_model(f) = Z_tissue_FEM(f) + Z_CPE_dl(f; dl_k, dl_alpha)

where ``Z_tissue_FEM`` comes from a single pure-Dirichlet FEM run
(``Results_NoInterface_Lempka_Spectrum``) shared by both scenarios.

Inputs (all over the same 200-point frequency grid):

- ``Results_LempkaImpedance_spectrum/`` — full FEM, Day 1 CPE
  (``dl_k=1.5e6, dl_alpha=0.8``).
- ``Results_LempkaImpedance_day15_spectrum/`` — full FEM, Day 15 CPE
  (``dl_k=1e6, dl_alpha=0.75``).
- ``Results_NoInterface_Lempka_Spectrum/`` — pure Dirichlet, shared
  tissue response.

Usage
-----
    ossdbs input_interface_lempka2009_spectrum.json
    ossdbs input_interface_lempka2009_day15_spectrum.json
    ossdbs input_no_interface_lempka2009_spectrum.json
    python plot_impedance_spectrum_day1_vs_day15.py
"""

import impedancefitter as ifit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_Z(path: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path)
    freq = df["freq"].to_numpy()
    Z = df["real"].to_numpy() + 1j * df["imag"].to_numpy()
    return freq, Z


# --- CPE_dl (double-layer) equivalent circuit model ---
ecm_cpe = ifit.get_equivalent_circuit_model("CPE_dl")

# --- Load FEM results ---
freq_day1, sim_Z_day1 = _load_Z("Results_LempkaImpedance_spectrum/impedance.csv")
freq_day15, sim_Z_day15 = _load_Z(
    "Results_LempkaImpedance_day15_spectrum/impedance.csv"
)
freq_tissue, sim_Z_tissue = _load_Z("Results_NoInterface_Lempka_Spectrum/impedance.csv")

assert np.all(np.isclose(freq_day1, freq_day15))
assert np.all(np.isclose(freq_day1, freq_tissue))
frequencies = freq_day1
omega = 2.0 * np.pi * frequencies

# --- Lempka 2009 CPE_dl parameters (must match the JSON inputs) ---
scenarios = {
    "Day 1": {
        "params": {"dl_k": 1.5e6, "dl_alpha": 0.8},
        "sim_Z": sim_Z_day1,
        "color": "C0",
    },
    "Day 15": {
        "params": {"dl_k": 1.0e6, "dl_alpha": 0.75},
        "sim_Z": sim_Z_day15,
        "color": "C1",
    },
}

# --- Reconstructed impedance: FEM tissue + ECM CPE_dl ---
for s in scenarios.values():
    s["Z_model"] = sim_Z_tissue + ecm_cpe.eval(omega=omega, **s["params"])

# --- Bode plot: full FEM vs decomposition, Day 1 and Day 15 overlaid ---
fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))

for name, s in scenarios.items():
    label_full = f"{name} - full FEM"
    label_model = f"{name} - FEM tissue + CPE_dl"
    ax_mag.plot(frequencies, np.abs(s["sim_Z"]), color=s["color"], label=label_full)
    ax_mag.plot(
        frequencies,
        np.abs(s["Z_model"]),
        color=s["color"],
        ls="dashed",
        label=label_model,
    )
    ax_phase.plot(
        frequencies,
        np.angle(s["sim_Z"], deg=True),
        color=s["color"],
        label=label_full,
    )
    ax_phase.plot(
        frequencies,
        np.angle(s["Z_model"], deg=True),
        color=s["color"],
        ls="dashed",
        label=label_model,
    )

ax_mag.plot(
    frequencies, np.abs(sim_Z_tissue), color="k", ls=":", label="FEM tissue only"
)
ax_phase.plot(
    frequencies,
    np.angle(sim_Z_tissue, deg=True),
    color="k",
    ls=":",
    label="FEM tissue only",
)

ax_mag.set_ylabel(r"|Z| / $\Omega$")
ax_mag.set_xscale("log")
ax_mag.set_yscale("log")
ax_mag.legend(fontsize=8)
ax_mag.grid(True, which="both", ls=":", alpha=0.5)
ax_mag.set_title("Day 1 vs Day 15 — Full FEM vs FEM-tissue + CPE_dl")

ax_phase.set_xlabel("Frequency / Hz")
ax_phase.set_ylabel("Phase / deg")
ax_phase.set_xscale("log")
ax_phase.grid(True, which="both", ls=":", alpha=0.5)

fig.tight_layout()
fig.savefig("impedance_spectrum_day1_vs_day15_bode.pdf")
fig.savefig("impedance_spectrum_day1_vs_day15_bode.svg")
plt.show()

# --- Relative error per scenario: |Z_model - Z_sim| / |Z_sim| ---
fig2, ax_err = plt.subplots(figsize=(7, 4))
for name, s in scenarios.items():
    rel_err = np.abs(s["Z_model"] - s["sim_Z"]) / np.abs(s["sim_Z"])
    ax_err.plot(frequencies, rel_err * 100.0, color=s["color"], label=name)
ax_err.set_xlabel("Frequency / Hz")
ax_err.set_ylabel("Relative error / %")
ax_err.set_xscale("log")
ax_err.legend()
ax_err.grid(True, which="both", ls=":", alpha=0.5)
ax_err.set_title("Decomposition error: |FEM tissue + CPE_dl - full FEM| / |full FEM|")
fig2.tight_layout()
fig2.savefig("impedance_spectrum_day1_vs_day15_error.pdf")
fig2.savefig("impedance_spectrum_day1_vs_day15_error.svg")
plt.show()

# --- Nyquist plot ---
fig3, ax3 = plt.subplots(figsize=(7, 6))
for name, s in scenarios.items():
    ax3.plot(
        s["sim_Z"].real,
        -s["sim_Z"].imag,
        color=s["color"],
        label=f"{name} - full FEM",
    )
    ax3.plot(
        s["Z_model"].real,
        -s["Z_model"].imag,
        color=s["color"],
        ls="dashed",
        label=f"{name} - FEM tissue + CPE_dl",
    )
ax3.set_xlabel(r"Re(Z) / $\Omega$")
ax3.set_ylabel(r"$-$Im(Z) / $\Omega$")
ax3.set_aspect("equal")
ax3.legend(fontsize=8)
ax3.grid(True, ls=":", alpha=0.5)
ax3.set_title("Nyquist — Day 1 vs Day 15")
fig3.tight_layout()
fig3.savefig("impedance_spectrum_day1_vs_day15_nyquist.pdf")
fig3.savefig("impedance_spectrum_day1_vs_day15_nyquist.svg")
plt.show()
