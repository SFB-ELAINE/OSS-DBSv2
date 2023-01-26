
from abc import ABC, abstractmethod


class Preconditioner(ABC):

    @abstractmethod
    def to_dictionary(self) -> dict:
        pass


class BDDCPreconditioner(Preconditioner):

    def __init__(self) -> None:
        self.type = 'bddc'
        self.coarsetype = 'local'

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
