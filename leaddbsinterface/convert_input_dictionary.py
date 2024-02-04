import argparse
import json

from ossdbs.utils.settings import Settings

from .default_settings import (
    initialize_default_settings,
    load_default_for_lead,
    update_default_dict,
)
from .lead_settings import LeadSettings


def main():
    """Lead-DBS/OSS-DBS Interface."""
    parser = argparse.ArgumentParser(
        prog="Lead-DBS to OSS-DBS converter",
        description="""Converts mat-files written by Lead-DBS
                                                    to OSS-DBS input dictionaries.""",
    )
    parser.add_argument(
        "leaddbs_dictionary",
        type=str,
        help="input dictionary in mat format provided by Lead-DBS",
    )
    parser.add_argument(
        "--hemi_side",
        type=int,
        choices=[0, 1],
        required=True,
        help="hemisphere side, 0 is right hemisphere, 1 is left hemisphere",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="",
        help="specify where to store the outputs of the OSS-DBS run",
    )
    args = parser.parse_args()

    # get default settings (alternatively, set using GUI)
    settings = initialize_default_settings()
    settings = load_default_for_lead(settings)

    # get settings from oss-dbs_parameters.mat
    ls = LeadSettings(args.leaddbs_dictionary)
    settings_leaddbs = ls.make_oss_settings(
        hemis_idx=args.hemi_side, output_path=args.output_path
    )

    # merge with default, overwrite if necessary
    settings = update_default_dict(settings, settings_leaddbs)

    # complete the input dictionary
    partial_settings = Settings(settings)
    complete_settings = partial_settings.complete_settings()

    # replace ending in input dictionaries
    new_input = args.leaddbs_dictionary.replace(".mat", ".json")

    # export json dict
    json_settings = json.dumps(complete_settings, indent=2)
    with open(new_input, "w") as outfile:
        outfile.write(json_settings)


if __name__ == "__main__":
    main()
