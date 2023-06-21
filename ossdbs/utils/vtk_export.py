from dataclasses import dataclass
import ngsolve


@dataclass
class FieldSolution:

    solution: ngsolve.CoefficientFunction
    label: str
    mesh: ngsolve.comp.Mesh
    complex: bool

    def save(self, filename: str) -> None:
        names = ["{}_real".format(self.label)]
        if complex:
            names.append("{}_imag".format(self.label))

        coefficients = [self.solution.real]
        if complex:
            coefficients.append(self.solution.imag)

        ngsolve.VTKOutput(ma=self.mesh,
                          coefs=coefficients,
                          names=names,
                          filename=filename,
                          subdivision=0
                          ).Do()
