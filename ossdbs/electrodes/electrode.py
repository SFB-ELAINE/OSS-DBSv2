from abc import ABC, abstractmethod
import netgen


class Electrode(ABC):
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
                 position: tuple = (0, 0, 0),
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

    def bounding_box(self):
        start, end = self.generate_geometry().bounding_box
        return tuple(start), tuple(end)
