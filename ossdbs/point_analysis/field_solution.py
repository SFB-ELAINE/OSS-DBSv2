
from dataclasses import dataclass
from ossdbs.fem import Solution
import ngsolve


@dataclass
class FieldSolution:

    solution: Solution
    mesh: ngsolve.comp.Mesh

    def save(self, filename: str) -> None:
        names = ["potential_real",
                 "potential_imag",
                 "current_density_real",
                 "current_density_imag"]

        coefficients = [self.solution.potential.real,
                        self.solution.potential.imag,
                        self.solution.current_density.real,
                        self.solution.current_density.imag]

        ngsolve.VTKOutput(ma=self.mesh,
                          coefs=coefficients,
                          names=names,
                          filename=filename,
                          subdivision=0
                          ).Do()
