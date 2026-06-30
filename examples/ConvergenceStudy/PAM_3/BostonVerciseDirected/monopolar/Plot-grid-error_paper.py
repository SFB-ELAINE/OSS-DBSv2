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
    "fine_hp_material_refinement": 16,
}

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
for pathway, label in zip(pathways_to_plot, pathway_labels, strict=False):
    pathway_label_dict[pathway] = label


convergence_threshold = 5.0  # in %

data = pd.read_csv("pam_results_summary.csv")

# drop the neuron mesh-size strategies and renumber to the canonical map
data = data[data["study_name"].isin(PAPER_ROMAN)].copy()
data["_order"] = data["study_name"].map(PAPER_ROMAN)
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

g = sns.PairGrid(data, x_vars=columns_to_plot, y_vars=["roman"], height=4)

# white = converged, red = above the 5% threshold (evaluated per error column)
palette = {True: "#FF0000", False: "#FFFFFF"}
for i, ax in enumerate(g.axes.flat):
    col = columns_to_plot[i]
    is_not_converged = (
        (data[col] > convergence_threshold) if i > 1 else [False] * len(data)
    )
    sns.stripplot(
        data=data,
        x=col,
        y="roman",
        hue=is_not_converged,
        palette=palette,
        ax=ax,
        size=10,
        orient="h",
        jitter=False,
        linewidth=1,
        edgecolor="black",
        legend=False,
    )
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#444444")
    ax.set(xlabel=labels[i], xscale=scales[i], ylabel="Strategy")
    # set consistent limit except for M1 upper extrem
    if i > 1:
        if i != 3:
            ax.set(xlim=(-0.1, 2.5))
        else:
            ax.set(xlim=(-0.1, None))
        ax.axvline(convergence_threshold, color="red", linestyle="--", alpha=0.5)
sns.despine(left=True, bottom=False)
plt.savefig("pam_convergence_overview_errors_paper.pdf")
plt.savefig("pam_convergence_overview_errors_paper.svg")
plt.close()
