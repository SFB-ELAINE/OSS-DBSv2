
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

    def append(self, contact: Contact) -> None:
        if contact.active:
            self.__active.append(contact)
        if contact.floating:
            self.__floating.append(contact)

    def set_voltage(self, contact_id: str, value: float) -> None:
        """Set voltage value of a named contact.

        Parameters
        ----------
        contact_id : str
            Name of the contact.

        value : float
            Voltage value.
        """
        for contact in self.__active:
            if contact.name is contact_id:
                contact.voltage = value

    def active(self) -> List[Contact]:
        """Returns a collection of all active contacts.

        Returns
        -------
        list of Contacts
        """
        return self.__active

    def floating(self) -> List[Contact]:
        """Returns a collection of all floating contacts.

        Returns
        -------
        list of Contacts
        """
        return self.__floating

    def current_values(self) -> dict:
        """Returns the current values of each active contact.

        Returns
        -------
        dict
        """
        return {contact.name: contact.current for contact in self.__active}

    def voltage_values(self) -> dict:
        """Returns the voltage values of each active contact.

        Returns
        -------
        dict
        """
        return {contact.name: contact.voltage for contact in self.__active}

    def floating_impedance_values(self) -> dict:
        """Returns the floating impedance values of each floating contact.

        Returns
        -------
        dict
        """
        return {contact.name: contact.surface_impedance
                for contact in self.__floating}

    def copy(self) -> 'Contacts':
        """Returns a copy of this collection.

        Returns
        -------
        Contacts
        """
        return Contacts(contacts=self.__active + self.__floating)
