"""Panel B: vary the boundary condition on the active contact.

Non-active segmented contacts are kept insulating (default Neumann BC)
and the case (BrainSurface) stays grounded at 0 V. Only C2's BC changes:

  B1_dirichlet  pure Dirichlet, V = 1 V, no surface impedance
  B2_mild_Z     V = 1 V + resistive surface impedance R = 1 kOhm
  B3_strong_Z   V = 1 V + resistive surface impedance R = 10 kOhm

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
MILD_R = 1e3
STRONG_R = 1e4


def make_contact_list(c2_scheme: str) -> list[dict]:
    """Build the contact list for one Panel B scheme."""
    c2 = {
        "Contact_ID": 2,
        "Active": True,
        "Current[A]": 0.0,
        "Voltage[V]": 1.0,
        "Floating": False,
    }
    if c2_scheme == "mild_Z":
        c2["SurfaceImpedance"] = {"Model": "R", "Parameters": {"R": MILD_R}}
    elif c2_scheme == "strong_Z":
        c2["SurfaceImpedance"] = {"Model": "R", "Parameters": {"R": STRONG_R}}
    elif c2_scheme == "dirichlet":
        pass
    else:
        raise ValueError(f"Unknown c2 scheme: {c2_scheme}")

    contacts = [c2]
    for cid in NON_ACTIVE_IDS:
        contacts.append(
            {
                "Contact_ID": cid,
                "Active": False,
                "Current[A]": 0.0,
                "Voltage[V]": 0.0,
                "Floating": False,
            }
        )
    return contacts


RUNS = [
    ("B1_dirichlet", "dirichlet"),
    ("B2_mild_Z", "mild_Z"),
    ("B3_strong_Z", "strong_Z"),
]


def main() -> None:
    """Run Panel B (B1_dirichlet, B2_mild_Z, B3_strong_Z) sequentially."""
    ossdbs.set_logger()
    logger = logging.getLogger("ossdbs")

    with open("base_settings.json") as fp:
        base = json.load(fp)

    for run_id, scheme in RUNS:
        d = copy.deepcopy(base)
        d["Electrodes"][0]["Contacts"] = make_contact_list(scheme)
        d["Surfaces"] = [{"Name": "BrainSurface", "Active": True, "Voltage[V]": 0.0}]
        d["OutputPath"] = f"Results_{run_id}"
        logger.info("=== Panel B run: %s (scheme=%s) ===", run_id, scheme)
        main_run(d)
        remove_file_handler(logger)


if __name__ == "__main__":
    main()
