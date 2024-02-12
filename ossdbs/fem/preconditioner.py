# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from abc import ABC, abstractmethod

_logger = logging.getLogger(__name__)


class Preconditioner(ABC):
    """Define a preconditioner."""

    @abstractmethod
    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        pass


class BDDCPreconditioner(Preconditioner):
    """BDDC preconditioner, highly efficient for higher-order FEM."""

    def __init__(self, coarsetype: str = "h1amg") -> None:
        _logger.debug(f"Initialized BDDC preconditioner with coarsetype: {coarsetype}")
        self.type = "bddc"
        self.coarsetype = coarsetype

    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        return {"type": self.type, "coarsetype": self.coarsetype}


class LocalPreconditioner(Preconditioner):
    """Jacobi preconditioner."""

    def __init__(self) -> None:
        _logger.debug("Initialized Jacobi preconditioner")
        self.type = "local"

    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        return {"type": self.type}


class MultigridPreconditioner(Preconditioner):
    """Geometric multigrid preconditioner."""

    def __init__(self) -> None:
        _logger.debug("Initialized multigrid preconditioner")
        self.type = "multigrid"

    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        return {"type": self.type}


class AMGPreconditioner(Preconditioner):
    """Algebraic multigrid preconditioner."""

    def __init__(self) -> None:
        _logger.debug("Initialized AMG preconditioner")
        self.type = "h1amg"

    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        return {"type": self.type}


class DirectPreconditioner(Preconditioner):
    """Direct solver."""

    def __init__(self, inverse: str = "") -> None:
        _logger.debug("Use direct solver as preconditioner")
        self.type = "direct"
        self.inverse = inverse

    def to_dictionary(self) -> dict:
        """Prepare dictionary to be passed to NGSolve."""
        if self.inverse == "":
            return {"type": self.type}
        return {"type": self.type, "inverse": self.inverse}
