# Copyright 2023, 2024 Konstantin Butenko, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from .colecole4 import ColeCole4Model, ColeColeParameters

WhiteMatterColeCole3Default = ColeColeParameters(
    alpha=np.array([0.1, 0.1, 0.3]),
    eps_delta=np.array([32.0, 100.0, 4.0e4]),
    eps_inf=4.0,
    sigma=0.265,
    tau=np.array([7.958e-12, 7.958e-9, 53.052e-6]),
)

GrayMatterColeCole3Default = ColeColeParameters(
    alpha=np.array([0.1, 0.15, 0.22]),
    eps_delta=np.array([45.0, 400.0, 2.0e5]),
    eps_inf=4.0,
    sigma=0.239,
    tau=np.array([7.958e-12, 15.915e-9, 106.103e-6]),
)

CSFColeCole3Default = ColeColeParameters(
    alpha=np.array([0.1, 0.0, 0.0]),
    eps_delta=np.array([65.0, 40.0, 0.0]),
    eps_inf=4.0,
    sigma=2.0,
    tau=np.array([7.96e-12, 1.592e-9, 159.155e-6]),
)


BloodColeCole3Default = ColeColeParameters(
    alpha=np.array([0.1, 0.1, 0.0]),
    eps_delta=np.array([56.0, 5200.0, 0.0]),
    eps_inf=4.0,
    sigma=0.7,
    tau=np.array([8.38e-12, 132.63e-9, 0]),
)

default_cole_cole3_parameters = {
    "Gray matter": GrayMatterColeCole3Default,
    "Unknown": GrayMatterColeCole3Default,
    "White matter": WhiteMatterColeCole3Default,
    "CSF": CSFColeCole3Default,
    "Blood": BloodColeCole3Default,
}


class ColeCole3Model(ColeCole4Model):
    """Cole-Cole model with three dispersions.

    Notes
    -----
    The model values are chosen as described
    in [Zimmermann2021]_.

    References
    ----------
    .. [Zimmermann2021] Zimmermann, J. and van Rienen, U. (2021)
                        Ambiguity in the interpretation of the low-frequency
                        dielectric properties of biological tissues.
                        Bioelectrochemistry, 140, 107773.
                        dx.doi.org/10.1016/j.bioelechem.2021.107773

    """

    def __init__(self, parameters: ColeColeParameters):
        self._parameters = parameters
        if not self._parameters.assert_order(3):
            raise ValueError(
                "ColeCole3Model requires information about three dispersions"
            )
