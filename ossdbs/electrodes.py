
from typing import List
from ossdbs.contacts import Contacts, Contact
from ossdbs.electrode_models.electrode import ElectrodeModel
from ossdbs.encapsulatin_layer import EncapsulatingLayers
import netgen


class Electrodes:

    def __init__(self,
                 electrodes: List[ElectrodeModel],
                 contacts: List[Contact]) -> None:

        self.__electrodes = electrodes
        self.__contacts = contacts

    def geometry(self) -> netgen.libngpy._NgOCC.TopoDS_Shape:
        geometry = sum([elec.geometry() for elec in self.__electrodes])
        names = [contact.name for contact in self.__contacts
                 if contact.active or contact.floating]
        for edge in geometry.edges:
            if edge.name in names:
                edge.maxh = 0.1
        return geometry

    def encapsulating_layer(self, thickness: float) -> EncapsulatingLayers:
        return EncapsulatingLayers(self.__electrodes, thickness)

    def contacts(self) -> List[Contact]:
        return Contacts([contact for contact in self.__contacts
                         if contact.active or contact.floating])
