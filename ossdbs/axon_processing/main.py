# Copyright 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
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
        "pathway_file", type=str, help="input file with pathway information"
    )
    parser.add_argument(
        "path_to_time_domain_solution",
        type=str,
        help="path where time domain solution is stored",
    )

    args = parser.parse_args()

    settings = {}
    settings["PathwayFile"] = args.pathway_file
    settings["OutputPath"] = args.path_to_time_domain_solution

    set_logger(level=args.loglevel)

    run_PAM(settings)

    _logger.info("Loading settings from input file")
