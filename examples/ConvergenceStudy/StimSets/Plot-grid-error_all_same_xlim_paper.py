import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- 1. Style Settings (Black Background / White Text) ---

plt.rcParams["text.usetex"] = True
plt.style.use("dark_background")  # Sets background black and basic elements white
plt.rcParams.update(
    {
        "axes.facecolor": "black",
        "figure.facecolor": "black",
        "axes.edgecolor": "white",
        "axes.labelcolor": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "grid.color": "gray",  # Subtle grid
        "text.color": "white",
        "axes.labelsize": 18,
        "text.latex.preamble": (
            r"\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother"  # noqa: E501
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
}

# StimSets paper: restrict to the strategies that were also run for the
# StimSets VTA study (minus the neuron mesh-size one). Canonical order.
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

pathway_label_dict = dict(zip(pathways_to_plot, pathway_labels, strict=False))
convergence_threshold = 5.0  # in %

# Load data
data = pd.read_csv("pam_results_summary.csv")

# keep the paper subset in canonical order and relabel the romans
data = data[data["study_name"].isin(STIMSETS_STRATEGIES)].copy()
data["_order"] = data["study_name"].map(STIMSETS_STRATEGIES.index)
data = data.sort_values("_order").drop(columns="_order")
data["roman"] = (
    data["study_name"].map(lambda s: rf"\rom{{{PAPER_ROMAN[s]}}}").astype("string")
)

# Build column lists
columns_to_plot = ["time", "dofs"]
labels = ["Time / s", "DOFs"]
scales = ["log", "log"]

for pathway in pathways_to_plot:
    pw_id = f"Pathway_status_{pathway}.json_activated_max_error"
    if pw_id not in data.columns:
        raise ValueError(f"Pathway {pathway} not in dataset.")
    columns_to_plot.append(pw_id)
    labels.append(rf"Err. {pathway_label_dict[pathway]} / %")
    scales.append("linear")

# --- 2. Plotting ---
# We do not pass 'hue' here because 'hue' will be column-specific
g = sns.PairGrid(data, x_vars=columns_to_plot, y_vars=["roman"], height=4)

# Define our highlight palette
# True (Not Converged) -> Red, False (Converged) -> White
palette = {True: "#FF0000", False: "#FFFFFF"}

# Iterate through axes and columns to apply specific logic
for i, ax in enumerate(g.axes.flat):
    col_name = columns_to_plot[i]

    # Logic: If it's an error column, calculate convergence per point
    # If it's time/dofs, we can default to False (White)
    if "max_error" in col_name:
        is_not_converged = data[col_name] > convergence_threshold
    else:
        is_not_converged = [False] * len(data)

    sns.stripplot(
        data=data,
        x=col_name,
        y="roman",
        hue=is_not_converged,
        palette=palette,
        ax=ax,
        size=10,
        orient="h",
        jitter=False,
        linewidth=1,
        edgecolor="black",
        legend=False,  # Remove individual legends
    )

    # --- 3. Formatting each axis ---
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#444444")
    ax.set(xlabel=labels[i], xscale=scales[i], ylabel="Strategy")

    if i > 1:  # Pathway error columns
        ax.set(xlim=(-0.3, 8))
        # Optional: Add a vertical line for the threshold
        ax.axvline(convergence_threshold, color="red", linestyle="--", alpha=0.5)

sns.despine(left=True, bottom=False)

plt.savefig("pam_convergence_overview_errors_all_same_paper.pdf")
plt.savefig("pam_convergence_overview_errors_all_same_paper.svg")
plt.close()
