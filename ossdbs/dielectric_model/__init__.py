# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Dielectric properties of required tissues."""
from .colecole3 import ColeCole3Model, default_cole_cole3_parameters
from .colecole4 import ColeCole4Model, ColeColeParameters, default_cole_cole4_parameters
from .constant import ConstantModel, ConstantParameters, default_constant_parameters
from .dielectric_model import DielectricModel

dielectric_models = {
    "ColeCole4": ColeCole4Model,
    "ColeCole3": ColeCole3Model,
    "Constant": ConstantModel,
}

dielectric_model_parameters = {
    "ColeCole4": ColeColeParameters,
    "ColeCole3": ColeColeParameters,
    "Constant": ConstantParameters,
}

default_dielectric_parameters = {
    "ColeCole4": default_cole_cole4_parameters,
    "ColeCole3": default_cole_cole3_parameters,
    "Constant": default_constant_parameters,
}


__all__ = (
    "DielectricModel",
    "ColeCole4Model",
    "ColeCole3Model",
    "ColeColeParameters",
    "ConstantModel",
    "ConstantParameters",
    "default_cole_cole4_parameters",
    "default_cole_cole3_parameters",
    "default_constant_parameters",
    "dielectric_models",
    "dielectric_model_parameters",
    "default_dielectric_parameters",
)
