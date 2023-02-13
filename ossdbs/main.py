
import ngsolve
import sys
import json
from ossdbs.impedance_analysis.impedance_analysis import impedance_analysis
from ossdbs.point_analysis.point_analysis import point_analysis


def main() -> None:

    path = sys.argv[1]
    with open(path, 'r') as json_file:
        input = json.load(json_file)

    if 'impedance' in sys.argv:
        with ngsolve.TaskManager():
            impedance_analysis(input)
        return

    with ngsolve.TaskManager():
        point_analysis(input)


if __name__ == '__main__':
    main()
