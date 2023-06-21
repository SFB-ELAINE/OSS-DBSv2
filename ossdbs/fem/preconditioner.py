
from abc import ABC, abstractmethod


class Preconditioner(ABC):

    @abstractmethod
    def to_dictionary(self) -> dict:
        pass


class BDDCPreconditioner(Preconditioner):

    def __init__(self, coarsetype: str = "h1amg") -> None:
        self.type = 'bddc'
        self.coarsetype = coarsetype

    def to_dictionary(self) -> dict:
        return {'type': self.type, 'coarsetype': self.coarsetype}


class LocalPreconditioner(Preconditioner):
    def __init__(self) -> None:
        self.type = 'local'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class MultigridPreconditioner(Preconditioner):
    def __init__(self) -> None:
        self.type = 'multigrid'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class AMGPreconditioner(Preconditioner):
    def __init__(self) -> None:
        self.type = 'h1amg'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class DirectPreconditioner(Preconditioner):
    def __init__(self, inverse: str = "") -> None:
        self.type = 'direct'
        self.inverse = inverse

    def to_dictionary(self) -> dict:
        if self.inverse == "":
            return {'type': self.type}
        return {'type': self.type, 'inverse': self.inverse}
