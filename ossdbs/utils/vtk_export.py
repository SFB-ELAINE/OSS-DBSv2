from dataclasses import dataclass

import ngsolve


@dataclass
class FieldSolution:
    solution: ngsolve.CoefficientFunction
    label: str
    mesh: ngsolve.comp.Mesh
    is_complex: bool

    def save(self, filename: str) -> None:
        names = [f"{self.label}_real"]
        if self.is_complex:
            names.append(f"{self.label}_imag")

        coefficients = [self.solution.real]
        if self.is_complex:
            coefficients.append(self.solution.imag)

        ngsolve.VTKOutput(
            ma=self.mesh,
            coefs=coefficients,
            names=names,
            filename=filename,
            subdivision=0,
        ).Do()
