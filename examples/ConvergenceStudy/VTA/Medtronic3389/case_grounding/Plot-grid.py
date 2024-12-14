import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Thanks to StackOverflow: https://stackoverflow.com/questions/34177378/pyplot-annotation-roman-numerals
# Turn on LaTeX formatting for text
plt.rcParams["text.usetex"] = True

# Place the command in the text.latex.preamble using rcParams
# ruff: noqa: E501
plt.rcParams["text.latex.preamble"] = (
    r"\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother"
)

data = pd.read_csv("vta_results_summary.csv")
data["roman"] = data["roman"].astype("string")

columns_to_plot = ["time", "dofs", "imp_rel_error", "ngs_vta_volume_rel_error"]
labels = ["Time / s", "DOFs", "Rel. err. impedance / %%", "Rel. err. VTA / %%"]
scales = ["linear", "log", "linear", "linear"]

g = sns.PairGrid(
    data, x_vars=data[columns_to_plot], y_vars=["roman"], height=10, aspect=0.2
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
# Use the same x axis limits on all columns and add better labels
# g.set(xlim=(0, 25), xlabel="Crashes", ylabel="")

for ax, label, scale in zip(g.axes.flat, labels, scales):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    # set labels and scales
    ax.set(xlabel=label)
    ax.set(xscale=scale)
    ax.set(ylabel="Strategy")
sns.despine(left=True, bottom=False)
plt.savefig("test.png")
plt.savefig("test.svg")
plt.close()
