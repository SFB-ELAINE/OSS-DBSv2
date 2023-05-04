import sys
import json
from ossdbs.impedance_analysis import impedance_analysis
from ossdbs.point_analysis import point_analysis
from ossdbs import set_logger


def main() -> None:

    path = sys.argv[1]
    with open(path, 'r') as json_file:
        input = json.load(json_file)

    if 'impedance' in sys.argv:
        impedance_analysis(input)
        return

    point_analysis(input)


if __name__ == '__main__':
    # default logger
    set_logger()
    main()
