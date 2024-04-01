import json
import logging

import ossdbs
import ossdbs.axon_processing

ossdbs.set_logger(level=logging.INFO)

pathway_file = "Allocated_axons_parameters.json"
time_domain_solution = "oss_time_result_PAM.h5"
with open(pathway_file) as fp:
    pathways_dict = json.load(fp)
print(pathways_dict)
model_type = pathways_dict["Axon_Model_Type"]

if "MRG2002" in model_type:
    downsampled = model_type == "MRG2002_DS"
    test = ossdbs.axon_processing.MRG2002(pathways_dict, "pam_results")
    print("Model is downsampled: ", downsampled)
    test.downsampled = downsampled
elif "McNeal1976" in model_type:
    test = ossdbs.axon_processing.McNeal1976(pathways_dict, "pam_results")
else:
    raise NotImplementedError(f"Model {model_type} not yet implemented.")

test.load_solution(time_domain_solution)

test.process_pathways(scaling=1.0, scaling_index=None)
