import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pathway = "DRTT_v1"
cutoff = 100

best_pathway_dir = "Results_PAM_best"
compare_pathway_dir = "Results_PAM_default"
time_result_file = "oss_time_result_PAM.h5"
activation_result_file = f"Axon_state_{pathway}.csv"

time_file_best = h5py.File(os.path.join(best_pathway_dir, time_result_file))
time_file_compare = h5py.File(os.path.join(compare_pathway_dir, time_result_file))
n_axons = len(time_file_best[pathway]["Status"])

activation_file_best = os.path.join(best_pathway_dir, activation_result_file)
activation_best = pd.read_csv(activation_file_best)
activation_file_compare = os.path.join(compare_pathway_dir, activation_result_file)
activation_compare = pd.read_csv(activation_file_compare)

if not os.path.isdir("figures"):
    os.mkdir("figures")

current_index = 0
for i in range(n_axons):
    n_neurons = np.count_nonzero(np.isclose(activation_compare["idx"], i))
    potential_best = np.array(time_file_best[pathway][f"axon{i}"]["Potential[V]"])
    potential_compare = np.array(time_file_compare[pathway][f"axon{i}"]["Potential[V]"])
    for idx_neuron in range(n_neurons):
        status_best = activation_best["status"].iloc[current_index]
        status_compare = activation_compare["status"].iloc[current_index]
        same_status = np.all(np.isclose(status_best, status_compare))
        if not same_status:
            p_best = potential_best[idx_neuron]
            p_compare = potential_compare[idx_neuron]
            plt.plot(p_best[:cutoff])
            plt.plot(p_compare[:cutoff])
            plt.savefig(f"figures/plot_{i}_{idx_neuron}.pdf")
            plt.close()
        current_index += 1

time_file_best.close()
time_file_compare.close()
