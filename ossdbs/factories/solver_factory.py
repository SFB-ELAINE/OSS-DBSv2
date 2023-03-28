
from ossdbs.fem import BDDCPreconditioner
from ossdbs.fem import LocalPreconditioner
from ossdbs.fem import MultigridPreconditioner
from ossdbs.fem import Solver, CGSolver, GMRESSolver


class SolverFactory:

    @classmethod
    def create(cls, parameters: dict) -> Solver:
        solver_type = parameters['Type']
        solver = {'CG': CGSolver, 'GMRES': GMRESSolver}[solver_type]
        preconditioner = {'bddc': BDDCPreconditioner(),
                          'local': LocalPreconditioner(),
                          'multigrid': MultigridPreconditioner()
                          }[parameters['Preconditioner']]

        return solver(precond_par=preconditioner,
                      printrates=parameters['PrintRates'],
                      maxsteps=parameters['MaximumSteps'],
                      precision=parameters['Precision'])
