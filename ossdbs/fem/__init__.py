from .mesh import Mesh
from .preconditioner import (LocalPreconditioner,
                             DirectPreconditioner,
                             BDDCPreconditioner,
                             AMGPreconditioner,
                             MultigridPreconditioner)
from .solver import Solver, GMRESSolver, CGSolver
from .volume_conductor import (ConductivityCF,
                               VolumeConductor,
                               VolumeConductorNonFloating,
                               VolumeConductorFloating,
                               VolumeConductorFloatingImpedance)

SOLVERS = {'CG': CGSolver,
           'GMRES': GMRESSolver}


PRECONDITIONERS = {'bddc': BDDCPreconditioner,
                   'local': LocalPreconditioner,
                   'multigrid': MultigridPreconditioner,
                   'h1amg': AMGPreconditioner,
                   'direct': DirectPreconditioner
                   }


__all__ = ['LocalPreconditioner',
           'DirectPreconditioner',
           'BDDCPreconditioner',
           'MultigridPreconditioner',
           'AMGPreconditioner',
           'Solver',
           'ConductivityCF',
           'GMRESSolver',
           'CGSolver',
           'Mesh',
           'VolumeConductor',
           'VolumeConductorNonFloating',
           'VolumeConductorFloating',
           'VolumeConductorFloatingImpedance']
