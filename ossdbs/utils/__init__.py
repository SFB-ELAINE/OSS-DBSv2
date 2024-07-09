# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""OSS-DBS utility functions."""
from .materials import have_dielectric_properties_changed
from .settings import Settings
from .type_check import TypeChecker
from .vtk_export import FieldSolution

__all__ = (
    "Settings",
    "TypeChecker",
    "FieldSolution",
    "have_dielectric_properties_changed",
)
