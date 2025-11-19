# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""FEM part of OSS-DBS."""

from .mesh import Mesh
from .preconditioner import (
    AMGPreconditioner,
    BDDCPreconditioner,
    DirectPreconditioner,
    LocalPreconditioner,
    MultigridPreconditioner,
)
from .solver import CGSolver, DirectSolver, GMRESSolver, Solver
from .volume_conductor import (
    ConductivityCF,
    VolumeConductor,
    VolumeConductorFloating,
    VolumeConductorFloatingImpedance,
    VolumeConductorNonFloating,
)

SOLVERS = {"CG": CGSolver, "GMRES": GMRESSolver, "Direct": DirectSolver}


PRECONDITIONERS = {
    "bddc": BDDCPreconditioner,
    "local": LocalPreconditioner,
    "multigrid": MultigridPreconditioner,
    "h1amg": AMGPreconditioner,
    "direct": DirectPreconditioner,
}


__all__ = [
    "AMGPreconditioner",
    "BDDCPreconditioner",
    "CGSolver",
    "ConductivityCF",
    "DirectPreconditioner",
    "DirectSolver",
    "GMRESSolver",
    "LocalPreconditioner",
    "Mesh",
    "MultigridPreconditioner",
    "Solver",
    "VolumeConductor",
    "VolumeConductorFloating",
    "VolumeConductorFloatingImpedance",
    "VolumeConductorNonFloating",
]
