"""Plot the StimSets truncation test (dark background, LaTeX labels).

Three panels, mirroring the single-protocol figure in
``PAM_3/BostonVerciseDirected/TruncationTest/plot_truncation.py`` but with the
error split into the mean and the worst case over the current protocols:

1. mean absolute activation error per pathway across the protocol sweep
2. maximum absolute activation error per pathway (the decisive panel)
3. time steps per unit solution and total H5 size of the eight unit solutions
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("dark_background")

plt.rcParams["text.usetex"] = True
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["legend.fontsize"] = 9
plt.rcParams["text.latex.preamble"] = r"\usepackage{siunitx}"

# Pathways with relevant activation across the protocol sweep. Unlike the
# single-protocol figure, M1 face is included -- it is silent for the one
# monopolar protocol but reaches 36 % mean activation over the sweep. Medial
# lemniscus stays out: 0.04 % mean activation, so its error is not meaningful.
PATHWAYS = [
    ("M1_cf_face_right", "M1 face"),
    ("M1_cf_lowerex_right", "M1 lower extr."),
    ("M1_cf_upperex_right", "M1 upper extr."),
    ("R_M1_hdp_face_right", "HDP M1 face"),
    ("R_M1_hdp_lowerex_right", "HDP M1 lower extr."),
    ("R_M1_hdp_upperex_right", "HDP M1 upper extr."),
    ("cerebellothalamic_right", "Cerebellothalamic"),
    ("gpe2stn_ass_right", "Pallido-subthal. Assoc"),
    ("gpe2stn_sm_right", "Pallido-subthal. Motor"),
]

MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*", "<"]

data = pd.read_csv("truncation_results_summary.csv")
data["truncation_ratio"] = data["truncation_ratio"].astype(str)

# The untruncated run is the reference: its error is zero by construction and
# it is not a point on the truncation-ratio axis (it corresponds to no
# truncation at all). Including it would suggest the error rises from "none"
# to ratio 5, so the error panels show only the truncated runs. The cost
# panel keeps it, where it is the meaningful baseline.
errors = data[data["truncation_ratio"] != "none"].reset_index(drop=True)

err_labels = errors["truncation_ratio"].values
err_pos = np.arange(len(err_labels))
cost_labels = data["truncation_ratio"].values
cost_pos = np.arange(len(cost_labels))

fig, (ax1, ax2, ax3) = plt.subplots(
    1, 3, figsize=(15, 4.6), gridspec_kw={"wspace": 0.35}, facecolor="black"
)


def plot_errors(ax, suffix, ylabel):
    """Draw one error panel over the truncated runs."""
    for (pathway, label), marker in zip(PATHWAYS, MARKERS, strict=False):
        ax.plot(
            err_pos,
            errors[f"{pathway}{suffix}"].values,
            marker=marker,
            label=label,
            linewidth=1.5,
            markersize=6,
        )
    ax.set_xlabel("Truncation ratio")
    ax.set_ylabel(ylabel)
    ax.set_xticks(err_pos)
    ax.set_xticklabels(err_labels)
    ax.set_ylim(bottom=-1)
    ax.axhline(0, color="white", linewidth=0.5, linestyle="--", alpha=0.3)


plot_errors(ax1, "_error", r"Mean absolute activation error / \si{\percent}")
plot_errors(ax2, "_max_error", r"Max. absolute activation error / \si{\percent}")

# --- Right panel: time steps and total unit-solution size ---
color_ts = "white"
color_fs = "C1"
bar_width = 0.6

ax3.bar(cost_pos, data["timesteps"].values, width=bar_width, color=color_ts, alpha=0.6)
ax3.set_xlabel("Truncation ratio")
ax3.set_ylabel("Time steps", color=color_ts)
ax3.tick_params(axis="y", labelcolor=color_ts)
ax3.set_xticks(cost_pos)
ax3.set_xticklabels(cost_labels)

ax3b = ax3.twinx()
ax3b.plot(
    cost_pos,
    data["h5_size_total_MB"].values / 1e3,
    marker="o",
    color=color_fs,
    linewidth=1.5,
    markersize=6,
)
ax3b.set_ylabel(r"Unit solutions / \si{\giga\byte}", color=color_fs)
ax3b.tick_params(axis="y", labelcolor=color_fs)

ax3b.spines["right"].set_color(color_fs)
ax3.spines["left"].set_color(color_ts)

# One shared legend below the panels: inside ax2 it would sit on top of the
# ratio-5 data, which reaches 100 %.
handles, labels = ax1.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="lower center",
    ncol=5,
    framealpha=0.2,
    bbox_to_anchor=(0.5, -0.13),
)

plt.savefig("truncation_test.pdf", bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig("truncation_test.svg", bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print("Saved truncation_test.pdf / .svg")
