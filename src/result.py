from src.output import Output
from dataclasses import dataclass
import ngsolve


@dataclass
class Result:

    mesh: ngsolve.Mesh
    potential: ngsolve.GridFunction

    def save(self, path: str = '') -> None:
        Output(mesh=self.mesh,
               potential=self.potential.real,
               field=-ngsolve.grad(self.potential).real * 1e3
               ).save(path)
