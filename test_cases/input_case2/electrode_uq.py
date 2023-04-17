import sys
import json
import ngsolve
import pandas as pd
from ossdbs.point_analysis import point_analysis

def uq_electrode():

    path = sys.argv[1]
    with open(path, 'r') as json_file:
       input = json.load(json_file)


    df = pd.read_csv("./SNEX-100_TC_Electrodes_geometry_wo_offset.csv")

    for row in range(len(df)):

        electrode_name = df._get_value(row, 'ELECTRODE_NAME')
        core_length = float(df._get_value(row, 'CORE_ELECTRODE_LENGTH')) * 1e-3
        core_diameter = float(df._get_value(row, 'CORE_ELECTRODE_DIAMETER')) * 1e-3
        core_tubing_length = float(df._get_value(row, 'CORE_TUBING_LENGTH')) * 1e-3
        core_tubing_diameter = float(df._get_value(row, 'CORE_TUBING_DIAMETER')) * 1e-3
        outer_electrode_length = float(df._get_value(row, 'OUTER_ELECTRODE_LENGTH')) * 1e-3
        outer_electrode_diameter = float(df._get_value(row, 'OUTER_ELECTRODE_DIAMETER')) * 1e-3
        outer_tubing_diameter = float(df._get_value(row, 'OUTER_TUBING_DIAMETER')) * 1e-3

        write_json(core_length,
                   core_diameter,
                   core_tubing_length,
                   core_tubing_diameter,
                   outer_electrode_length,
                   outer_electrode_diameter,
                    outer_tubing_diameter  )

        input["OutputPath"] = "case2_" + str(row)

        with ngsolve.TaskManager():
            point_analysis(input)


def write_json(core_length,
               core_diameter,
               core_tubing_length,
               core_tubing_diameter,
               outer_electrode_length,
               outer_electrode_diameter,
               outer_tubing_diameter  ):
    data = {
        "CoreElectrodeLength[mm]": core_length,
        "CoreElectrodeDiameter[mm]": core_diameter,
        "CoreTubingLength[mm]": core_tubing_length,
        "CoreTubingDiameter[mm]": core_tubing_diameter,
        "OuterElectrodeLength[mm]": outer_electrode_length,
        "OuterElectrodeDiameter[mm]": outer_electrode_diameter,
        "OuterTubingDiameter[mm]": outer_tubing_diameter,
        "LeadLength[mm]": 100.0
    }

    with open('../ossdbs/electrodes/electrode_models/custom/micro_probes_SNEX100_custom.json', 'w') as file:
        json.dump(data, file)


if __name__ == '__main__':
    uq_electrode()