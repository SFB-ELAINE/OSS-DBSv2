import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Thanks to StackOverflow: https://stackoverflow.com/questions/34177378/pyplot-annotation-roman-numerals
# Turn on LaTeX formatting for text
plt.rcParams["text.usetex"] = True
plt.rcParams["axes.labelsize"] = 18

# Place the command in the text.latex.preamble using rcParams
# ruff: noqa: E501
plt.rcParams["text.latex.preamble"] = (
    r"\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother"
)

# Worst-case companion to Plot-grid.py: instead of the mean over the 254
# protocols, look at the worst protocol per strategy (max rel. volume error
# and min Dice). Thresholds are looser than the mean plot because a single
# protocol's error is always larger than the average.
max_error_threshold = 2.0  # max. rel. VTA volume error, in %
min_dice_threshold = 0.98  # min. Dice vs. best

data = pd.read_csv("vta_results_summary.csv")
data["roman"] = data["roman"].astype("string")

# The summary already carries the reference run as a row labelled "Best"
# (rel. error 0, Dice 1). Relative errors are stored as fractions, so
# convert to percent for plotting.
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
    linewidth=1,
    edgecolor="w",
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
    ax.yaxis.grid(True)
    # set labels and scales
    ax.set(xlabel=label)
    ax.set(xscale=scale)
    ax.set(ylabel="Strategy")
sns.despine(left=True, bottom=False)
plt.savefig("vta_convergence_overview_max.pdf")
plt.savefig("vta_convergence_overview_max.svg")
plt.close()
