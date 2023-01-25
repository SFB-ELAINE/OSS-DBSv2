
from abc import ABC, abstractmethod


class PreconditionerParameters(ABC):

    @abstractmethod
    def to_dictionary(self) -> dict:
        pass


class BDDCPreconditioner(PreconditionerParameters):

    def __init__(self) -> None:
        self.type = 'bddc'
        self.coarsetype = 'local'

    def to_dictionary(self) -> dict:
        return {'type': self.type, 'coarsetype': self.coarsetype}


class LocalPreconditioner(PreconditionerParameters):
    def __init__(self) -> None:
        self.type = 'local'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class MultigridPreconditioner(PreconditionerParameters):
    def __init__(self) -> None:
        self.type = 'multigrid'

    def to_dictionary(self) -> dict:
        return {'type': self.type}
