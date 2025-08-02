# Copyright 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import json
import logging

from ossdbs import set_logger
from ossdbs.api import run_PAM

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
        "input_dictionary",
        type=str,
        help="input dictionary in JSON format as used for preparation of StimSets run",
    )
    parser.add_argument(
        "--scaling", type=float, help="specify scalar scaling of solution", default=1.0
    )
    parser.add_argument(
        "--scaling_index",
        type=int,
        help="specify index of the scaled solution",
        default=None,
    )

    args = parser.parse_args()
    set_logger(level=args.loglevel)

    _logger.info("Loading settings from input file")
    _logger.debug(f"Input file: {args.input_dictionary}")
    with open(args.input_dictionary) as json_file:
        input_settings = json.load(json_file)

    if not input_settings["StimSets"]["Active"]:
        _logger.info("No StimSets will be run.")
    input_settings["Scaling"] = args.scaling
    input_settings["ScalingIndex"] = args.scaling_index
    input_settings["CurrentVector"] = None

    run_PAM(input_settings)
