import json
import logging
import multiprocessing as mp
from pathlib import Path

import ossdbs
import ossdbs.axon_processing


def main() -> None:
    """Run the bundled McNeal1976 PAM example."""
    ossdbs.set_logger(level=logging.INFO)

    base_dir = Path(__file__).resolve().parent.parent
    pathway_file = Path(__file__).resolve().parent / "Allocated_axons_parameters.json"
    time_domain_solution = base_dir / "Results_rh" / "oss_time_result_PAM.h5"
    with open(pathway_file) as fp:
        pathways_dict = json.load(fp)
    print(pathways_dict)
    model_type = pathways_dict["Axon_Model_Type"]

    if "MRG2002" in model_type:
        downsampled = model_type == "MRG2002_DS"
        test = ossdbs.axon_processing.MRG2002(
            pathways_dict, str(Path(__file__).resolve().parent / "pam_results")
        )
        print("Model is downsampled: ", downsampled)
        test.downsampled = downsampled
    elif "McNeal1976" in model_type:
        test = ossdbs.axon_processing.McNeal1976(
            pathways_dict, str(Path(__file__).resolve().parent / "pam_results")
        )
    else:
        raise NotImplementedError(f"Model {model_type} not yet implemented.")

    td_solution = test.load_solution(str(time_domain_solution))
    test.process_pathways(td_solution, scaling=1.0, scaling_index=None)


if __name__ == "__main__":
    mp.freeze_support()
    main()
