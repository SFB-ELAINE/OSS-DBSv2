import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- 1. Style Settings (Black Background / White Text) ---
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
        # Explicitly uniform font sizes
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    }
)

# Canonical strategy -> Roman numeral for the paper figures. Shared verbatim
# with the StimSets VTA figure (Plot-grid-max_paper.py) so a numeral always
# denotes the same mesh-refinement strategy across both figures.
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

# StimSets paper: restrict to the strategies also run for the StimSets VTA
# study, so the PAM and VTA figures show the exact same strategy set.
STIMSETS_STRATEGIES = [
    "default",
    "fine",
    "material_refinement",
    "edge_single_material_refinement",
    "hp_refinement",
    "fine_hp_material_refinement",
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
    "CST M1 face",
    "CST M1 lower extr.",
    "CST M1 upper extr.",
    "HDP M1 face",
    "HDP M1 lower extr.",
    "HDP M1 upper extr.",
    "Dentatorubrothalamic",
    "GPe -> STN assoc.",
    "GPe -> STN motor",
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
    data["study_name"].map(lambda s: ROMAN_NUMERALS[PAPER_ROMAN[s]]).astype("string")
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
    labels.append(rf"{pathway_label_dict[pathway]} / %")
    scales.append("linear")

# --- 2. Plotting ---
# We do not pass 'hue' here because 'hue' will be column-specific
g = sns.PairGrid(data, x_vars=columns_to_plot, y_vars=["roman"], height=4, aspect=1)

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
        ax.set(xlim=(0, 17))
        # Optional: Add a vertical line for the threshold
        ax.axvline(convergence_threshold, color="red", linestyle="--", alpha=0.5)
    else:  # Time / DOFs columns
        ax.set(xscale="log")
        if i == 0:
            ax.set(xlim=(10, 1000))
        else:
            ax.set(xlim=(1e4, 2e6))

sns.despine(left=True, bottom=False)

plt.savefig("pam_convergence_overview_errors_all_same_paper.pdf")
plt.savefig("pam_convergence_overview_errors_all_same_paper.svg")
plt.close()
