"""Plot truncation test results.

Two-panel figure:
  Left:  per-pathway absolute activation error vs truncation ratio
  Right: number of time steps vs truncation ratio
Includes the untruncated reference as a data point.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["text.usetex"] = True
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["legend.fontsize"] = 10
plt.rcParams["text.latex.preamble"] = r"\usepackage{siunitx}"

# Pathways with nonzero activation (skip M1_face, medial_lemniscus)
pathways = [
    ("M1_cf_lowerex_right", "M1 lower extr."),
    ("M1_cf_upperex_right", "M1 upper extr."),
    ("R_M1_hdp_face_right", "HDP M1 face"),
    ("R_M1_hdp_lowerex_right", "HDP M1 lower extr."),
    ("R_M1_hdp_upperex_right", "HDP M1 upper extr."),
    ("cerebellothalamic_right", "Cerebellothalamic"),
    ("gpe2stn_ass_right", "Pallido-subthal. Assoc"),
    ("gpe2stn_sm_right", "Pallido-subthal. Motor"),
]

data = pd.read_csv("truncation_results_summary.csv")
data["truncation_ratio"] = data["truncation_ratio"].astype(str)

# Evenly-spaced x positions so the "none" point isn't squashed at the left
x_labels = data["truncation_ratio"].values
x_pos = np.arange(len(x_labels))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"wspace": 0.35})

# --- Left panel: activation error ---
markers = ["o", "s", "^", "v", "D", "P", "X", "*"]
for (pw, label), marker in zip(pathways, markers, strict=False):
    errors = data[f"{pw}_error"].values
    ax1.plot(x_pos, errors, marker=marker, label=label, linewidth=1.5, markersize=6)

ax1.set_xlabel("Truncation ratio")
ax1.set_ylabel(r"Absolute activation error / \si{\percent}")
ax1.set_xticks(x_pos)
ax1.set_xticklabels(x_labels)
ax1.legend(loc="upper right", framealpha=0.9)
ax1.set_ylim(bottom=-1)
ax1.axhline(0, color="grey", linewidth=0.5, linestyle="--")

# --- Right panel: time steps and file size ---
color_ts = "C0"
color_fs = "C3"
bar_width = 0.6

ax2.bar(x_pos, data["timesteps"].values, width=bar_width, color=color_ts, alpha=0.8)
ax2.set_xlabel("Truncation ratio")
ax2.set_ylabel("Time steps", color=color_ts)
ax2.tick_params(axis="y", labelcolor=color_ts)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(x_labels)

# Twin axis for file size
ax2b = ax2.twinx()
ax2b.plot(
    x_pos,
    data["h5_size_MB"].values,
    marker="o",
    color=color_fs,
    linewidth=1.5,
    markersize=6,
)
ax2b.set_ylabel(r"H5 file size / \si{\mega\byte}", color=color_fs)
ax2b.tick_params(axis="y", labelcolor=color_fs)

plt.savefig("truncation_test.pdf", bbox_inches="tight")
plt.savefig("truncation_test.svg", bbox_inches="tight")
plt.close()
print("Saved truncation_test.pdf / .svg")
