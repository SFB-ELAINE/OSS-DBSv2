from abc import ABC, abstractmethod
import netgen


class AbstractElectrode(ABC):
    """Deep Brain Simulation electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    translation : tuple
        Translation vector (x,y,z) of electrode.f

    Methods
    -------
    generate_geometry()
        Generate mesh of electrode.
    """
    @abstractmethod
    def __init__(elf,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0)) -> None:
        pass

    @abstractmethod
    def generate_geometry() -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        pass
