

from dataclasses import dataclass
from typing import List
from ossdbs.electrodes.electrode import Electrode
from ossdbs.encapsulatin_layer import EncapsulatingLayers


@dataclass
class ElectrodeContact:
    name: str = ''
    active: bool = False
    current: float = 0.0
    floating: bool = False
    voltage: float = 0.0
    surface_impedance: float = 0.0


class Electrodes:

    def __init__(self, electrodes: list, contacts: list) -> None:

        self.__electrodes = electrodes
        self.__contacts = contacts

    def geometry(self):
        geometry = sum([electrode.geometry()
                        for electrode in self.__electrodes])

        names = [contact.name for contact in self.__contacts
                 if contact.active or contact.floating]
        for edge in geometry.edges:
            if edge.name in names:
                edge.maxh = 0.1
        return geometry

    def encapsulating_layer(self, thickness: float) -> EncapsulatingLayers:
        return EncapsulatingLayers(self.__electrodes, thickness)

    def contacts(self) -> List[ElectrodeContact]:
        return [contact for contact in self.__contacts
                if contact.active or contact.floating]

    def active_contacts(self) -> List[str]:
        return sorted([contact.name for contact in self.__contacts
                       if contact.active])

    def floating_contacts(self) -> List[str]:
        return sorted([contact.name for contact in self.__contacts
                       if contact.floating and not contact.active])

    def current_values(self) -> dict:
        return {contact.name: contact.current
                for contact in self.__contacts if contact.active}

    def voltage_values(self) -> dict:
        return {contact.name: contact.voltage
                for contact in self.__contacts if contact.active}

    def floating_impedance_values(self) -> dict:
        return {contact.name: contact.surface_impedance
                for contact in self.__contacts
                if contact.floating and not contact.active}

    def capsule_geometries(self) -> List:
        return [electrode.capsule_geometry(thickness=0.1, max_h=0.1)
                for electrode in self.__electrodes]

    def electrode_geometries(self) -> List:
        return [electrode.geometry() for electrode in self.__electrodes]
