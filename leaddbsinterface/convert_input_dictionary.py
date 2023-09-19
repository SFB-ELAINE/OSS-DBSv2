import argparse
import json

from .default_settings import load_default_for_lead
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

    # get settings from oss-dbs_parameters.mat
    ls = LeadSettings(args.leaddbs_dictionary)
    settings = ls.make_oss_settings(
        hemis_idx=args.hemi_side, output_path=args.output_path
    )

    # add default settings (alternatively, set using GUI)
    settings = load_default_for_lead(settings)

    # replace ending in input dictionaries
    new_input = args.leaddbs_dictionary.replace(".mat", ".json")

    # export json dict
    json_settings = json.dumps(settings, indent=2)
    with open(new_input, "w") as outfile:
        outfile.write(json_settings)


if __name__ == "__main__":
    main()
