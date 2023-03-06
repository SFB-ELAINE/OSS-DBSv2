
from dataclasses import dataclass
from typing import List


@dataclass
class Contact:
    name: str = ''
    active: bool = False
    current: float = 0.0
    floating: bool = False
    voltage: float = 0.0
    surface_impedance: complex = 0.0j


class Contacts:

    def __init__(self, contacts: List[Contact]) -> None:
        self.__active = [contact for contact in contacts if contact.active]
        self.__floating = [contact for contact in contacts
                           if contact.floating and not contact.active]

    def set_voltage(self, contact_id: str, value: float) -> None:
        for contact in self.__active:
            if contact.name is contact_id:
                contact.voltage = value

    def active_contacts(self) -> List[Contact]:
        return self.__active

    def floating_contacts(self) -> List[Contact]:
        return self.__floating

    def current_values(self) -> dict:
        return {contact.name: contact.current for contact in self.__active}

    def voltage_values(self) -> dict:
        return {contact.name: contact.voltage for contact in self.__active}

    def floating_impedance_values(self) -> dict:
        return {contact.name: contact.surface_impedance
                for contact in self.__floating}

    def copy(self) -> 'Contacts':
        return Contacts(contacts=self.__active+self.__floating)
