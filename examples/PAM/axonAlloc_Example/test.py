import logging

import ossdbs
from ossdbs.axon_processing import AxonModels

ossdbs.set_logger(level=logging.INFO)

stim_dir = "./"
hemis_idx = 1
description_file = "oss-dbs_parameters.mat"

axon_model = AxonModels(stim_dir, hemis_idx, description_file)
axon_model.convert_fibers_to_axons()
