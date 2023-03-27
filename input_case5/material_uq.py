import sys
import json
import ngsolve
import pandas as pd
from ossdbs.point_analysis import point_analysis

def uq_material():

    path = sys.argv[1]
    with open(path, 'r') as json_file:
       input = json.load(json_file)


    df = pd.read_csv("./parameters_GM_CC4.csv")

    for row in range(len(df)):

        alpha_1 = df._get_value(row, 'Alpha_1')
        alpha_2 = df._get_value(row, 'Alpha_2')
        alpha_3 = df._get_value(row, 'Alpha_3')
        alpha_4 = df._get_value(row, 'Alpha_4')
        alpha = [alpha_1, alpha_2, alpha_3, alpha_4]

        epsilon_delta_1 = df._get_value(row, 'EpsilonDelta_1')
        epsilon_delta_2 = df._get_value(row, 'EpsilonDelta_2')
        epsilon_delta_3 = df._get_value(row, 'EpsilonDelta_3')
        epsilon_delta_4 = df._get_value(row, 'EpsilonDelta_4')
        epsilon_delta = [epsilon_delta_1, epsilon_delta_2, epsilon_delta_3, epsilon_delta_4]

        epsilon_inf = float(df._get_value(row, 'EpsilonInfinite'))

        sigma = float(df._get_value(row, 'Sigma'))

        tau_1 = df._get_value(row, 'Tau_1')
        tau_2 = df._get_value(row, 'Tau_2')
        tau_3 = df._get_value(row, 'Tau_3')
        tau_4 = df._get_value(row, 'Tau_4')
        tau = [tau_1, tau_2, tau_3, tau_4]

        write_json(alpha,
                   epsilon_delta,
                   epsilon_inf,
                   sigma,
                   tau)

        input["OutputPath"] = "case5_" + str(row)

        with ngsolve.TaskManager():
            point_analysis(input)


def write_json(alpha,
               epsilon_delta,
               epsilon_inf,
               sigma,
               tau):
    data = {
        "Alpha": alpha,
        "EpsilonDelta": epsilon_delta,
        "EpsilonInfinite": epsilon_inf,
        "Sigma": sigma,
        "Tau": tau
    }

    with open('../ossdbs/dielectric_model/custom/GrayMatterColeCole4Model.json', 'w') as file:
        json.dump(data, file)


if __name__ == '__main__':
    uq_material()