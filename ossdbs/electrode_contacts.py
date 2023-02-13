
from typing import List
from dataclasses import dataclass


@dataclass
class ElectrodeContact:
    name: str = ''
    active: bool = False
    current: float = 0.0
    floating: bool = False
    voltage: float = 0.0
    surface_impedance: float = 0.0


class ContactCollection:

    def __init__(self, contacts: list = None) -> None:
        self.__contacts = contacts if contacts else []

    def append(self, contact: ElectrodeContact) -> None:
        self.__contacts.append(contact)

    def floating(self) -> List[str]:
        return sorted([contact.name for contact in self.__contacts
                       if contact.floating and not contact.active])

    def set_inactive(self, name: str) -> None:
        for contact in self.__contacts:
            if contact.name == name:
                contact.active = False

    def active(self) -> List[str]:
        return sorted([contact.name for contact in self.__contacts
                       if contact.active])

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
