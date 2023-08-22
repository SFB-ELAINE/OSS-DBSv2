import sys
import json
import ngsolve
import pandas as pd
from ossdbs.point_analysis import point_analysis


def uq_material():

    df = pd.read_csv("./parameters_GM_CC4.csv")

    for row in range(len(df)):

        alpha_1 = float(df._get_value(row, 'Alpha_1'))
        alpha_2 = float(df._get_value(row, 'Alpha_2'))
        alpha_3 = float(df._get_value(row, 'Alpha_3'))
        alpha_4 = float(df._get_value(row, 'Alpha_4'))
        alpha = [alpha_1, alpha_2, alpha_3, alpha_4]

        epsilon_delta_1 = float(df._get_value(row, 'EpsilonDelta_1'))
        epsilon_delta_2 = float(df._get_value(row, 'EpsilonDelta_2'))
        epsilon_delta_3 = float(df._get_value(row, 'EpsilonDelta_3'))
        epsilon_delta_4 = float(df._get_value(row, 'EpsilonDelta_4'))
        epsilon_delta = [epsilon_delta_1, epsilon_delta_2, epsilon_delta_3, epsilon_delta_4]

        epsilon_inf = float(df._get_value(row, 'EpsilonInfinite'))

        sigma = float(df._get_value(row, 'Sigma'))

        tau_1 = float(df._get_value(row, 'Tau_1'))
        tau_2 = float(df._get_value(row, 'Tau_2'))
        tau_3 = float(df._get_value(row, 'Tau_3'))
        tau_4 = float(df._get_value(row, 'Tau_4'))
        tau = [tau_1, tau_2, tau_3, tau_4]

        write_json(alpha,
                   epsilon_delta,
                   epsilon_inf,
                   sigma,
                   tau)

        path = sys.argv[1]
        with open(path, 'r') as json_file:
            input = json.load(json_file)

        input["OutputPath"] = "case5_" + str(row)

        with ngsolve.TaskManager():
            point_analysis(input)


def write_json(alpha,
               epsilon_delta,
               epsilon_inf,
               sigma,
               tau):
    data = {
        "Gray matter": {
            "Alpha": alpha,
            "EpsilonDelta": epsilon_delta,
            "EpsilonInfinite": epsilon_inf,
            "Sigma": sigma,
            "Tau": tau
        },
        "White matter": {
            "Alpha": [0.1, 0.15, 0.22, 0.0],
            "EpsilonDelta": [45.0, 400.0, 2.0e5, 4.5e7],
            "EpsilonInfinite": 4.0,
            "Sigma": 0.02,
            "Tau": [7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3]
        },
        "CSF": {
            "Alpha": [0.1, 0.15, 0.22, 0.0],
            "EpsilonDelta": [45.0, 400.0, 2.0e5, 4.5e7],
            "EpsilonInfinite": 4.0,
            "Sigma": 0.02,
            "Tau": [7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3]
        },
        "Blood": {
            "Alpha": [0.1, 0.15, 0.22, 0.0],
            "EpsilonDelta": [45.0, 400.0, 2.0e5, 4.5e7],
            "EpsilonInfinite": 4.0,
            "Sigma": 0.02,
            "Tau": [7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3]
        }
    }

    with open('./ColeCole4ModelCustom.json', 'w') as file:
        json.dump(data, file)


if __name__ == '__main__':
    uq_material()
