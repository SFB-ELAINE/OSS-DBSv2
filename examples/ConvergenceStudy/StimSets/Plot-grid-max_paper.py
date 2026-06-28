import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- Style Settings (Black Background / White Text) ---
# ruff: noqa: E501
plt.rcParams["text.usetex"] = True
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
        "axes.labelsize": 18,
        "text.latex.preamble": (
            r"\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother"
        ),
    }
)

# Canonical strategy -> Roman numeral for the paper figures. The neuron
# mesh-size strategies (default_meshsize, edge_meshsize) are excluded.
# Shared verbatim across the PAM and VTA convergence figures so a numeral
# always denotes the same mesh-refinement strategy everywhere.
PAPER_ROMAN = {
    "default": 1,
    "fine": 2,
    "very_fine": 3,
    "edge_refinement": 4,
    "fine_edge_refinement": 5,
    "very_fine_edge_refinement": 6,
    "edge_voxel_refinement": 7,
    "material_refinement": 8,
    "edge_single_material_refinement": 9,
    "edge_double_material_refinement": 10,
    "fine_edge_single_material_refinement": 11,
    "fine_edge_double_material_refinement": 12,
    "hp_refinement": 13,
    "hp_material_refinement": 14,
    "hp_double_material_refinement": 15,
}

# StimSets paper: restrict to the strategies that were also run for the
# StimSets PAM study (minus the neuron mesh-size one), so the StimSets PAM
# and VTA figures show the exact same strategy set. Canonical order.
STIMSETS_STRATEGIES = [
    "default",
    "fine",
    "fine_edge_refinement",
    "edge_single_material_refinement",
    "hp_refinement",
    "hp_material_refinement",
    "hp_double_material_refinement",
]

# Worst-case companion to Plot-grid_paper.py: the worst protocol per
# strategy (max rel. volume error and min Dice) instead of the mean.
max_error_threshold = 2.0  # max. rel. VTA volume error, in %
min_dice_threshold = 0.98  # min. Dice vs. best

data = pd.read_csv("vta_results_summary.csv")

# keep the paper subset + the reference row, in canonical order, and
# relabel the romans consistently (reference -> "Best")
order = [*STIMSETS_STRATEGIES, "best"]
data = data[data["strategy"].isin(order)].copy()
data["_order"] = data["strategy"].map(order.index)
data = data.sort_values("_order").drop(columns="_order")
data["roman"] = (
    data["strategy"]
    .map(lambda s: "Best" if s == "best" else rf"\rom{{{PAPER_ROMAN[s]}}}")
    .astype("string")
)

# Relative errors are stored as fractions, so convert to percent.
data["max_err_pct"] = data["max_rel_volume_error"] * 100.0

columns_to_plot = [
    "time",
    "DOF",
    "mean_VTA_volume_mm3",
    "max_err_pct",
    "min_dice_vs_best",
]
labels = [
    "Time / s",
    "DOFs",
    r"Mean VTA volume / mm$^3$",
    r"Max. rel. err. VTA / \%",
    "Min. Dice",
]
scales = ["log", "log", "linear", "linear", "linear"]

converged = (data["max_err_pct"] <= max_error_threshold) & (
    data["min_dice_vs_best"] >= min_dice_threshold
)

g = sns.PairGrid(data, x_vars=data[columns_to_plot], y_vars=["roman"], height=4)
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
