# Copyright 2023, 2024 Konstantin Butenko, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from ossdbs.electrodes.defaults import default_electrode_parameters


def imp_coord(settings):
    """Calculates the lead implantation coordinate for a given settings dict."""
    electrode = settings["Electrodes"][settings["ModelSide"]]
    tip_pos = np.array(
        [
            electrode["TipPosition"]["x[mm]"],
            electrode["TipPosition"]["y[mm]"],
            electrode["TipPosition"]["z[mm]"],
        ]
    )

    # Find offset between lead implantation coordinate at center of contact
    # and the tip position of the electrode
    elec_params = default_electrode_parameters[electrode["Name"]]
    direction = np.array(
        [
            electrode["Direction"]["x[mm]"],
            electrode["Direction"]["y[mm]"],
            electrode["Direction"]["z[mm]"],
        ]
    )
    unit_direction = direction / np.linalg.norm(direction)
    return tip_pos + elec_params.offset * unit_direction
