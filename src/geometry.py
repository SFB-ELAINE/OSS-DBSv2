
from abc import ABC, abstractmethod


class Geometry(ABC):

    @abstractmethod
    def generate_mesh(self):
        pass
