import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

stim_sets_combinations = 254
pathways = [
    "M1_cf_face_right",
    "M1_cf_lowerex_right",
    "M1_cf_upperex_right",
    "R_M1_hdp_face_right",
    "R_M1_hdp_lowerex_right",
    "R_M1_hdp_upperex_right",
    "cerebellothalamic_right",
    "medial_lemniscus_right",
    "gpe2stn_ass_right",
    "gpe2stn_sm_right",
]

result_directories = [
    "Results_PAM_default",
    "Results_PAM_fine",
    "Results_PAM_fine_edge_refinement",
    "Results_PAM_material_refinement",
    "Results_PAM_edge_single_material_refinement",
]


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


best_result = pd.read_csv("Results_PAM_best_overview.csv")
columns = best_result.columns
print(type(columns), type(columns[0]))
for pathway in pathways:
    print(pathway in columns)

# get best result
mean_activations = []
mean_diff_dict = {}
max_diff_dict = {}
mean_diff_dict["roman"] = []
max_diff_dict["roman"] = []
for pw_idx, _ in enumerate(pathways):
    mean_diff_dict["roman"].append(rf"\rom{pw_idx + 1}")
    max_diff_dict["roman"].append(rf"\rom{pw_idx + 1}")
columns_to_plot = []
for result_directory in result_directories:
    print("++++++++++++++")
    print(result_directory)
    print("--------------")
    print("Pathway, Mean activation, Mean difference, Max. difference")
    df = pd.read_csv(result_directory + "_overview.csv")
    mean_diff_dict[result_directory] = []
    max_diff_dict[result_directory] = []
    columns_to_plot.append(result_directory)
    for pathway in pathways:
        best_activation = best_result[pathway].to_numpy()
        current_activation = df[pathway].to_numpy()
        difference = np.abs(best_activation - current_activation)
        mean_activation = np.mean(current_activation)
        std_activation = np.std(current_activation)
        mean_diff = np.mean(difference)
        std_diff = np.std(difference)
        max_diff = np.max(difference)
        print(f"{pathway}, {mean_activation:.2f}, {mean_diff:.2f}, {max_diff:.2f}")
        mean_diff_dict[result_directory].append(mean_diff)
        max_diff_dict[result_directory].append(max_diff)
    print("++++++++++++++")


labels = ["I", "II", "III", "IV", "V"]

for idx, data_dict in enumerate([mean_diff_dict, max_diff_dict]):
    x_max = 0
    data = pd.DataFrame(data_dict)
    for result_directory in result_directories:
        data["not_converged"] = data[result_directory] > convergence_threshold
        print(data[[result_directory, "not_converged"]])
        x_max = max(x_max, data[result_directory].max())

    g = sns.PairGrid(
        data,
        x_vars=data[columns_to_plot],
        y_vars=["roman"],
        height=4,
        hue="not_converged",
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

    for ax, _ in zip(g.axes.flat, labels):
        # Make the grid horizontal instead of vertical
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)
        # set labels and scales
        # ax.set(xlabel=label)
        ax.set(xlim=(-0.25, x_max * 1.1))
        ax.set(xlabel="Rel. error")
        ax.set(ylabel="Pathway")
    sns.despine(left=True, bottom=False)
    plt.savefig(f"pam_convergence_overview_{idx}.pdf")
    plt.savefig(f"pam_convergence_overview_{idx}.svg")
    plt.close()
