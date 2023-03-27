
from typing import List
from .contacts import Contacts
from .electrode_models.electrode import ElectrodeModel
from .encapsulatin_layer import EncapsulatingLayers
import netgen


class Electrodes:
    """Collection of electrode models.

    Attributes
    ----------
    electrodes : list of ElectrodeModel
        Collection of electrode models.

    contacts : list of Contact
        Collection of electrode contacts.
    """

    def __init__(self,
                 electrodes: List[ElectrodeModel],
                 contacts: Contacts) -> None:
        self.__electrodes = electrodes
        self.__contacts = contacts
        self.__max_h = 0.1

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        """Generate Geometry of all electrodes.

        Returns
        -------
        netgen.libngpy._NgOCC.TopoDS_Shape
        """

        geometry = sum([elec.geometry() for elec in self.__electrodes])
        names = [contact.name for contact
                 in self.__contacts.active() + self.__contacts.floating()]
        for edge in geometry.edges:
            if edge.name in names:
                edge.maxh = self.__max_h
        return geometry

    def encapsulating_layer(self, thickness: float) -> EncapsulatingLayers:
        """Return Collection of electrode encapsulating layers.

        Returns
        -------
        EncapsulatingLayers
        """
        return EncapsulatingLayers(self.__electrodes, thickness)

    def contacts(self) -> Contacts:
        """Return collection of active contacts and contacts of property
        floating

        Returns
        -------
        Contacts
        """
        return self.__contacts

    def set_max_h(self, max_h: float = 0.1) -> None:
        """Set the maximum height for a mesh element in these geometry of
        electrodes.

        Parameters
        ----------
        max_h : float
            Maximum height for a mesh element.
        """
        self.__max_h = max_h
