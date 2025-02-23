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

pathways_to_plot = [
    "M1_cf_lowerex_right",
    "M1_cf_upperex_right",
    # "M1_cf_face_right",  # not activated
    "R_M1_hdp_lowerex_right",
    "R_M1_hdp_upperex_right",
    "R_M1_hdp_face_right",
    "gpe2stn_sm_right",
    "gpe2stn_ass_right",
    "cerebellothalamic_right",
    # "medial_lemniscus_right"  # not activated
]

pathway_labels = [
    "M1 lower extr.",
    "M1 upper extr.",
    # "M1 face",  # not activated
    "HDP M1 lower extr.",
    "HDP M1 upper extr.",
    "HDP M1 face",
    "Pallido-subthalamic Motor",
    "Pallido-subthalamic Assoc",
    "Cerebellothalamic",
    # "Medial lemniscus"  # not activated
]

pathway_label_dict = {}
for pathway, label in zip(pathways_to_plot, pathway_labels):
    pathway_label_dict[pathway] = label


data = pd.read_csv("pam_results_summary.csv")
data["roman"] = data["roman"].astype("string")

for pathway in pathways_to_plot:
    pw_id = f"Pathway_status_{pathway}.json_activated"
    if pw_id not in data.columns:
        raise ValueError(f"Pathway {pathway} not in dataset.")

columns_to_plot = ["time", "dofs"]
labels = ["Time / s", "DOFs"]
scales = ["linear", "log"]
data["not_converged"] = False
for pathway in pathways_to_plot:
    pw_id = f"Pathway_status_{pathway}.json_activated"
    columns_to_plot.append(pw_id)
    labels.append(f"{pathway_label_dict[pathway]} / %")
    scales.append("linear")


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
