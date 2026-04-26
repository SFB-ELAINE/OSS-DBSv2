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

convergence_threshold = 3.0  # in %

data = pd.read_csv("vta_results_summary.csv")
data["roman"] = data["roman"].astype("string")

# Add "Best" row for absolute values
best_df = data.filter(regex="best", axis=1).iloc[0].to_dict()
best_df["time"] = best_df["best_time"]
best_df["dofs"] = best_df["best_dofs"]
best_df["ngs_vta_volume"] = best_df["best_ngs_vta_volume"]
best_df["imp"] = best_df["best_impedance"]
best_df["imp_rel_error"] = 0.0
best_df["ngs_vta_volume_rel_error"] = 0.0
best_df["roman"] = "Best"
best_df = pd.DataFrame([best_df])
data = pd.concat([data, best_df], ignore_index=True)
data.fillna(0.0, inplace=True)

columns_to_plot = [
    "time",
    "dofs",
    "ngs_vta_volume",
    "imp",
    "imp_rel_error",
    "ngs_vta_volume_rel_error",
]
labels = [
    "Time / s",
    "DOFs",
    r"VTA volume / mm$^3$",
    r"Impedance / $\Omega$",
    r"Rel. err. impedance / \%",
    r"Rel. err. VTA / \%",
]
scales = ["log", "log", "linear", "linear", "linear", "linear"]

converged = (data["imp_rel_error"] <= convergence_threshold) & (
    data["ngs_vta_volume_rel_error"] <= convergence_threshold
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
plt.savefig("vta_convergence_overview.pdf")
plt.savefig("vta_convergence_overview.svg")
plt.close()
