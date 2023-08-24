from abc import ABC, abstractmethod
import logging

_logger = logging.getLogger(__name__)


class Preconditioner(ABC):

    @abstractmethod
    def to_dictionary(self) -> dict:
        pass


class BDDCPreconditioner(Preconditioner):

    def __init__(self, coarsetype: str = "h1amg") -> None:
        _logger.debug("Initialized BDDC preconditioner with coarsetype: {}".format(coarsetype))
        self.type = 'bddc'
        self.coarsetype = coarsetype

    def to_dictionary(self) -> dict:
        return {'type': self.type, 'coarsetype': self.coarsetype}


class LocalPreconditioner(Preconditioner):
    def __init__(self) -> None:
        _logger.debug("Initialized Jacobi preconditioner")
        self.type = 'local'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class MultigridPreconditioner(Preconditioner):
    def __init__(self) -> None:
        _logger.debug("Initialized multigrid preconditioner")
        self.type = 'multigrid'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class AMGPreconditioner(Preconditioner):
    def __init__(self) -> None:
        _logger.debug("Initialized AMG preconditioner")
        self.type = 'h1amg'

    def to_dictionary(self) -> dict:
        return {'type': self.type}


class DirectPreconditioner(Preconditioner):
    def __init__(self, inverse: str = "") -> None:
        _logger.debug("Use direct solver as preconditioner")
        self.type = 'direct'
        self.inverse = inverse

    def to_dictionary(self) -> dict:
        if self.inverse == "":
            return {'type': self.type}
        return {'type': self.type, 'inverse': self.inverse}
