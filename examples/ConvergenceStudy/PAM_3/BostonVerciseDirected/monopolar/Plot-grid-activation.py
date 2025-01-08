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

convergence_threshold = 5.0  # in %

data = pd.read_csv("pam_results_summary.csv")
data["roman"] = data["roman"].astype("string")

columns_to_plot = ["time", "dofs"]
labels = ["Time / s", "DOFs"]
scales = ["linear", "log"]
data["not_converged"] = False
for column in data.columns:
    if "activated" in column:
        if "rel_err" in column or "best" in column:
            continue
        columns_to_plot.append(column)
        pathway_name = (
            column.replace("Pathway_status_", "")
            .replace(".json", "")
            .replace("_activated", "")
        )
        labels.append(rf"Activation {pathway_name} / \%")
        scales.append("linear")

        data["not_converged"] = (
            data["not_converged"] | data[column] > convergence_threshold
        )

g = sns.PairGrid(
    data, x_vars=data[columns_to_plot], y_vars=["roman"], height=4, hue="not_converged"
)
g.map(
    sns.stripplot,
    size=10,
    orient="h",
    jitter=False,
    palette="flare_r",
    linewidth=1,
    edgecolor="w",
)

for ax, label, scale in zip(g.axes.flat, labels, scales):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    # set labels and scales
    ax.set(xlabel=label)
    ax.set(xscale=scale)
    ax.set(ylabel="Strategy")
sns.despine(left=True, bottom=False)
plt.savefig("pam_convergence_overview.pdf")
plt.savefig("pam_convergence_overview.svg")
plt.close()
