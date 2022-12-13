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
    def __init__(self,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 translation: tuple = (0, 0, 0),
                 # contact_values: list = None
                 ) -> None:
        pass

    @abstractmethod
    def generate_geometry() -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        geometry : netgen.libngpy._NgOCC.TopoDS_Shape
        """
        pass

    @abstractmethod
    def rename_boundaries(self, boundary_names: dict) -> None:
        """Rename boundary names of electrode.

        Prameters
        ---------
        boundary_names : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        pass
