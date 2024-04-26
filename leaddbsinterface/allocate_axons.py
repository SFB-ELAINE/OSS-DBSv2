# Copyright 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import logging

from ossdbs import set_logger
from ossdbs.axon_processing import AxonModels


def main():
    """Axon allocation as specified by Lead-DBS."""
    parser = argparse.ArgumentParser(
        prog="Axon allocator based on Lead-DBS input",
        description="Converts mat-files written by Lead-DBS "
        "with OSS-DBS simulation results to "
        "an axon distribution used in pathway "
        "activation modeling.",
    )
    parser.add_argument(
        "stimulation_folder",
        type=str,
        help="full path to the stimulation folder",
    )
    parser.add_argument(
        "--hemi_side",
        type=int,
        choices=[0, 1],
        required=True,
        help="hemisphere side, 0 is right hemisphere, 1 is left hemisphere",
    )
    parser.add_argument(
        "--description_file",
        type=str,
        required=True,
        help="the file with Lead-DBS parameter information "
        "(usually oss-dbs_parameters.mat)",
    )

    parser.add_argument(
        "--loglevel", type=int, help="specify verbosity of logger", default=logging.INFO
    )

    args = parser.parse_args()
    set_logger(level=args.loglevel)

    # process Lead-DBS input
    axons_for_PAM = AxonModels(
        args.stimulation_folder, args.hemi_side, args.description_file
    )
    axons_for_PAM.convert_fibers_to_axons()


if __name__ == "__main__":
    main()
