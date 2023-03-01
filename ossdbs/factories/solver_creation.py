
from ossdbs.preconditioner import BDDCPreconditioner
from ossdbs.preconditioner import LocalPreconditioner
from ossdbs.preconditioner import MultigridPreconditioner
from ossdbs.solver import CGSolver, GMRESSolver


class SolverFactory:

    @classmethod
    def create(cls, parameters):
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
