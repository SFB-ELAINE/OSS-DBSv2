"""Panel A: vary the boundary condition on the non-active segmented contacts.

Active contact (C2 at +1 V) and case (BrainSurface, 0 V) are identical
across the three runs. The seven non-active segmented contacts
(C1, C3, C4, C5, C6, C7, C8) switch between:

  A1_floating    pure floating (no surface impedance)
  A2_insulating  natural Neumann BC (Active=false, Floating=false)
  A3_mild_Z      floating + resistive surface impedance R = 1 kOhm

Sequential execution per repository CLAUDE.md rule.
"""

import copy
import json
import logging

import ossdbs
from ossdbs.main import main_run


def remove_file_handler(logger):
    """Drop file handlers so each run logs to its own file."""
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)


NON_ACTIVE_IDS = [1, 3, 4, 5, 6, 7, 8]
MILD_R = 1e3  # Ohm, specific impedance — same scale as input_interface_1kOhm.json


def make_contact_list(scheme: str) -> list[dict]:
    """Build the contact list for one Panel A scheme."""
    contacts = [
        {
            "Contact_ID": 2,
            "Active": True,
            "Current[A]": 0.0,
            "Voltage[V]": 1.0,
            "Floating": False,
        }
    ]
    for cid in NON_ACTIVE_IDS:
        c = {
            "Contact_ID": cid,
            "Active": False,
            "Current[A]": 0.0,
            "Voltage[V]": 0.0,
            "Floating": False,
        }
        if scheme == "floating":
            c["Floating"] = True
        elif scheme == "mild_Z":
            c["Floating"] = True
            c["SurfaceImpedance"] = {"Model": "R", "Parameters": {"R": MILD_R}}
        elif scheme == "insulating":
            pass
        else:
            raise ValueError(f"Unknown scheme: {scheme}")
        contacts.append(c)
    return contacts


RUNS = [
    ("A1_floating", "floating"),
    ("A2_insulating", "insulating"),
    ("A3_mild_Z", "mild_Z"),
]


def main() -> None:
    """Run Panel A (A1_floating, A2_insulating, A3_mild_Z) sequentially."""
    ossdbs.set_logger()
    logger = logging.getLogger("ossdbs")

    with open("base_settings.json") as fp:
        base = json.load(fp)

    for run_id, scheme in RUNS:
        d = copy.deepcopy(base)
        d["Electrodes"][0]["Contacts"] = make_contact_list(scheme)
        d["Surfaces"] = [{"Name": "BrainSurface", "Active": True, "Voltage[V]": 0.0}]
        d["OutputPath"] = f"Results_{run_id}"
        # bddc fails to converge with floating contacts; matches the
        # solver block of examples/Electrochemistry/input_interface_1kOhm_floating.json
        if scheme in ("floating", "mild_Z"):
            d["Solver"]["Preconditioner"] = "local"
            d["Solver"]["MaximumSteps"] = 10000
        logger.info("=== Panel A run: %s (scheme=%s) ===", run_id, scheme)
        main_run(d)
        remove_file_handler(logger)


if __name__ == "__main__":
    main()
