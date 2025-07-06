# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Core module for volume-conductor model."""

from .conductivity import ConductivityCF
from .floating import VolumeConductorFloating
from .floating_impedance import VolumeConductorFloatingImpedance
from .nonfloating import VolumeConductorNonFloating
from .volume_conductor_model import VolumeConductor

__all__ = [
    "ConductivityCF",
    "VolumeConductor",
    "VolumeConductorFloating",
    "VolumeConductorFloatingImpedance",
    "VolumeConductorNonFloating",
]
