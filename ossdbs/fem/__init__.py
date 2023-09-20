from .mesh import Mesh
from .preconditioner import (
    AMGPreconditioner,
    BDDCPreconditioner,
    DirectPreconditioner,
    LocalPreconditioner,
    MultigridPreconditioner,
)
from .solver import CGSolver, GMRESSolver, DirectSolver, Solver
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
    "LocalPreconditioner",
    "DirectPreconditioner",
    "BDDCPreconditioner",
    "MultigridPreconditioner",
    "AMGPreconditioner",
    "Solver",
    "ConductivityCF",
    "GMRESSolver",
    "CGSolver",
    "Mesh",
    "VolumeConductor",
    "VolumeConductorNonFloating",
    "VolumeConductorFloating",
    "VolumeConductorFloatingImpedance",
]
