import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- Style Settings (Black Background / White Text) ---
plt.style.use("dark_background")  # black background, white foreground elements
plt.rcParams.update(
    {
        "axes.facecolor": "black",
        "figure.facecolor": "black",
        "axes.edgecolor": "white",
        "axes.labelcolor": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "grid.color": "gray",
        "text.color": "white",
        # Explicitly uniform font sizes (matching the PAM error figure)
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    }
)

# Canonical strategy -> Roman numeral for the paper figures. Shared verbatim
# with the StimSets PAM figure (Plot-grid-error_all_same_xlim_paper.py) so a
# numeral always denotes the same mesh-refinement strategy across both figures.
PAPER_ROMAN = {
    "default": 1,
    "fine": 2,
    "material_refinement": 3,
    "edge_single_material_refinement": 4,
    "hp_refinement": 5,
    "fine_hp_material_refinement": 6,
}

# A simple tuple lookup (0 is padded since Roman numerals start at 1)
ROMAN_NUMERALS = (
    "",
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
    "XVI",
)

# StimSets paper: restrict to the strategies also run for the StimSets PAM
# study, so the PAM and VTA figures show the exact same strategy set.
STIMSETS_STRATEGIES = [
    "default",
    "fine",
    "material_refinement",
    "edge_single_material_refinement",
    "hp_refinement",
    "fine_hp_material_refinement",
]

# Worst-case companion to Plot-grid_paper.py: the worst protocol per
# strategy (max rel. volume error and min Dice) instead of the mean.
max_error_threshold = 5.0  # max. rel. VTA volume error, in %
min_dice_threshold = 1.0  # min. Dice vs. best

data = pd.read_csv("vta_results_summary.csv")

# keep the paper subset + the reference row, in canonical order, and
# relabel the romans consistently (reference -> "Benchmark")
order = [*STIMSETS_STRATEGIES, "best"]
data = data[data["strategy"].isin(order)].copy()
data["_order"] = data["strategy"].map(order.index)
data = data.sort_values("_order").drop(columns="_order")
data["roman"] = (
    data["strategy"]
    .map(lambda s: "Benchmark" if s == "best" else ROMAN_NUMERALS[PAPER_ROMAN[s]])
    .astype("string")
)

# Relative errors are stored as fractions, so convert to percent.
data["max_err_pct"] = data["max_rel_volume_error"] * 100.0

columns_to_plot = [
    "time",
    "DOF",
    "max_err_pct",
    "min_dice_vs_best",
]
labels = [
    "Time / s",
    "DOFs",
    r"Max. relative error / %",
    "Min. Dice",
]
scales = ["log", "log", "linear", "linear"]

converged = (data["max_err_pct"] <= max_error_threshold) & (
    data["min_dice_vs_best"] >= min_dice_threshold
)

g = sns.PairGrid(data, x_vars=columns_to_plot, y_vars=["roman"], height=4, aspect=1)
g.map(
    sns.stripplot,
    size=10,
    orient="h",
    jitter=False,
    color="white",
    linewidth=1,
    edgecolor="black",
)

# Overlay orange dots for converged strategies
for ax_idx, col in enumerate(columns_to_plot):
    ax = g.axes.flat[ax_idx]
    for y_pos, (_, row) in enumerate(data.iterrows()):
        if converged.iloc[y_pos]:
            ax.scatter(
                row[col],
                y_pos,
                color="orange",
                s=100,
                zorder=5,
                edgecolor="w",
                linewidth=1,
            )

for ax, label, scale in zip(g.axes.flat, labels, scales, strict=False):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#444444")
    # set labels and scales
    ax.set(xlabel=label)
    ax.set(xscale=scale)
    ax.set(ylabel="Strategy")
sns.despine(left=True, bottom=False)
plt.savefig("vta_convergence_overview_max_paper.pdf")
plt.savefig("vta_convergence_overview_max_paper.svg")
plt.close()
