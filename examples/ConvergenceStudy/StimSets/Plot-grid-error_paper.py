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
}

# StimSets paper: restrict to the strategies that were also run for the
# StimSets VTA study (minus the neuron mesh-size one), so the StimSets PAM
# and VTA figures show the exact same strategy set. Canonical order.
STIMSETS_STRATEGIES = [
    "default",
    "fine",
    "fine_edge_refinement",
    "edge_single_material_refinement",
    "hp_refinement",
    "hp_material_refinement",
]

pathways_to_plot = [
    "M1_cf_face_right",
    "M1_cf_lowerex_right",
    "M1_cf_upperex_right",
    "R_M1_hdp_face_right",
    "R_M1_hdp_lowerex_right",
    "R_M1_hdp_upperex_right",
    "cerebellothalamic_right",
    "gpe2stn_ass_right",
    "gpe2stn_sm_right",
    "medial_lemniscus_right",
]

pathway_labels = [
    "M1 face",
    "M1 lower extr.",
    "M1 upper extr.",
    "HDP M1 face",
    "HDP M1 lower extr.",
    "HDP M1 upper extr.",
    "Cerebellothalamic",
    "Pallido-subthalamic Assoc",
    "Pallido-subthalamic Motor",
    "Medial lemniscus",
]

pathway_label_dict = {}
for pathway, label in zip(pathways_to_plot, pathway_labels, strict=False):
    pathway_label_dict[pathway] = label


convergence_threshold = 5.0  # in %

data = pd.read_csv("pam_results_summary.csv")

# keep the paper subset in canonical order and relabel the romans
data = data[data["study_name"].isin(STIMSETS_STRATEGIES)].copy()
data["_order"] = data["study_name"].map(STIMSETS_STRATEGIES.index)
data = data.sort_values("_order").drop(columns="_order")
data["roman"] = (
    data["study_name"].map(lambda s: rf"\rom{{{PAPER_ROMAN[s]}}}").astype("string")
)

for pathway in pathways_to_plot:
    pw_id = f"Pathway_status_{pathway}.json_activated_rel_error"
    if pw_id not in data.columns:
        raise ValueError(f"Pathway {pathway} not in dataset.")

columns_to_plot = ["time", "dofs"]
labels = ["Time / s", "DOFs"]
scales = ["log", "log"]
data["not_converged"] = False
for pathway in pathways_to_plot:
    pw_id = f"Pathway_status_{pathway}.json_activated_rel_error"
    columns_to_plot.append(pw_id)
    labels.append(rf"Err. {pathway_label_dict[pathway]} / \%")
    scales.append("linear")

    data["not_converged"] = data["not_converged"] | data[pw_id] > convergence_threshold

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

ax_count = 0
for ax, label, scale in zip(g.axes.flat, labels, scales, strict=False):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    # set labels and scales
    ax.set(xlabel=label)
    ax.set(xscale=scale)
    ax.set(ylabel="Strategy")
    # set consistent limit except for M1 face / cerebellothalamic
    if ax_count > 1:
        if ax_count not in (2, 8):  # M1 face, cerebellothalamic
            ax.set(xlim=(-0.1, 2.5))
        else:
            ax.set(xlim=(-0.1, None))
    ax_count += 1
sns.despine(left=True, bottom=False)
plt.savefig("pam_convergence_overview_errors_paper.pdf")
plt.savefig("pam_convergence_overview_errors_paper.svg")
plt.close()
