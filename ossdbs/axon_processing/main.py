# Copyright 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import json
import logging

from ossdbs import set_logger
from ossdbs.api import run_PAM
from ossdbs.utils.settings import Settings
from ossdbs.utils.type_check import TypeChecker

_logger = logging.getLogger(__name__)


def main() -> None:
    """Pathway activation function."""
    parser = argparse.ArgumentParser(
        prog="OSS-DBS pathway activation modelling.",
        description="Welcome to the pathway activation "
        "modelling submodule of OSS-DBS v2.",
        epilog="Please report bugs and errors on GitHub",
    )
    parser.add_argument(
        "--loglevel", type=int, help="specify verbosity of logger", default=logging.INFO
    )
    parser.add_argument(
        "input_dictionary", type=str, help="input dictionary in JSON format"
    )
    parser.add_argument(
        "pathway_file", type=str, help="input file with pathway information"
    )
    parser.add_argument(
        "--scaling", type=float, help="specify scalar scaling of solution", default=1.0
    )
    parser.add_argument(
        "--scaling_index", type=int, help="specify index of the scaled solution", default=None
    )
    parser.add_argument(
        "path_to_time_domain_solution",
        type=str,
        help="path where time domain solution is stored",
    )

    args = parser.parse_args()

    _logger.info("Loading settings from input file")
    _logger.debug(f"Input file: {args.input_dictionary}")
    with open(args.input_dictionary) as json_file:
        input_settings = json.load(json_file)

    settings = Settings(input_settings).complete_settings()
    TypeChecker.check(settings)

    settings["PathwayFile"] = args.pathway_file
    settings["OutputPath"] = args.path_to_time_domain_solution
    settings["Scaling"] = args.scaling
    settings["ScalingVector"] = None  # do not allow scaling vectors from terminal
    settings["ScalingIndex"] = args.scaling_index

    print(settings["ScalingVector"])

    set_logger(level=args.loglevel)

    run_PAM(settings)

    _logger.info("Loading settings from input file")
