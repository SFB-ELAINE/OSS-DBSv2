
from ossdbs.bounding_box import BoundingBox
from abc import ABC, abstractmethod
import netgen


class ElectrodeModel(ABC):
    """Deep Brain Simulation electrode.

    Attributes
    ----------
    rotation : float
        Rotation angle in degree of electrode.

    direction : tuple
        Direction vector (x,y,z) of electrode.

    position : tuple
        Position vector (x,y,z) of electrode tip.
    """
    @abstractmethod
    def __init__(self,
                 rotation: float = 0,
                 direction: tuple = (0, 0, 1),
                 position: tuple = (0, 0, 0),
                 ) -> None:
        pass

    @abstractmethod
    def geometry() -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of electrode.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """
        pass

    @abstractmethod
    def set_contact_names(self, contact_names: dict) -> None:
        """Set the names of electrode contacts.

        Parameters
        ----------
        contact_names : dict
            {'Body': 'body_name',
             'Contact_1': 'contact_name',
             'Contact_2': ...}
        """
        pass

    @abstractmethod
    def capsule_geometry(self, thickness: float, max_h: float) \
            -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate geometry of capsule layer around electrode.

        Parameters
        ----------
        thickness : float
            Thickness of encapsulating layer.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """

    def bounding_box(self) -> BoundingBox:
        """Return bounding box of electrode geometry.

        Returns
        -------
        BoundingBox
        """
        start, end = self.geometry().bounding_box
        return BoundingBox(tuple(start), tuple(end))
