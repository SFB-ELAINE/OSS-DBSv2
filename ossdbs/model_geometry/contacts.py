from dataclasses import dataclass
from typing import List
import logging
_logger = logging.getLogger(__name__)


@dataclass
class Contact:
    """Electrode contact settings.

    Notes
    -----
    This class stores the main parameters of the
    electrode contacts.

    General property:

        * `name`: Will be set during the geometry creation process.

    Mesh related properties:

        * `max_h`: Maximum element size on surface
        * `edge_max_h`: Maximum element size on contact edges

    Volume conductor model related properties:

        * `active`: Whether it needs a Dirichlet BC.
        * `floating`: Whether the voltage shall be fixed but unknown.
        * `current`: Assigned or computed current value.
        * `voltage`: Assigned or computed voltage value.
        * `surface_impedance`: Assigned surface impedance value.
    """

    name: str
    max_h: float = 1e10  # Netgen default
    edge_max_h: float = 1e10
    active: bool = False
    floating: bool = False
    current: float = 0.0
    voltage: float = 0.0
    surface_impedance: complex = 0.0j


def check_contact(contact: Contact):
    if contact.active and contact.floating:
        raise ValueError(
            "The contact {} has multiple roles. Please make sure that contacts are either active, floating or none of the two.".format(
                contact.name
            )
        )


class Contacts:
    """Wrapper class to classify contacts.

    Notes
    -----
    This class is intended to take the list of contacts of the model geometry and
    detect active, floating and unused contacts.
    """

    def __init__(self, contacts: List[Contact]) -> None:
        for contact in contacts:
            check_contact(contact)
        self._all_contacts = contacts
        self._active = [contact for contact in contacts if contact.active]
        self._floating = [contact for contact in contacts if contact.floating]
        self._unused = [
            contact
            for contact in contacts
            if not contact.floating and not contact.active
        ]

    def append(self, contact: Contact) -> None:
        self._all_contacts.append(contact)
        if contact.active:
            self._active.append(contact)
        elif contact.floating:
            self._floating.append(contact)
        else:
            self._unused.append(contact)

    @property
    def active(self) -> List[Contact]:
        """List of all active contacts.

        Returns
        -------
        list of Contacts
        """
        return self._active

    @property
    def floating(self) -> List[Contact]:
        """List of all floating contacts.

        Returns
        -------
        list of Contacts
        """
        return self._floating

    @property
    def unused(self) -> List[Contact]:
        """List of all unused (i.e. not active and not floating) contacts.

        Returns
        -------
        list of Contacts
        """
        return self._unused

    @property
    def currents(self) -> dict:
        """Returns the current values of each contact.

        Returns
        -------
        dict
        """
        return {contact.name: contact.current for contact in self._all_contacts}

    @currents.setter
    def currents(self, current_values: dict) -> None:
        """Set current values of contacts.

        Parameters
        ----------
        current_values : dict
            Current values. Not all contacts have to be present in the
            dictionary.

        """
        _logger.debug(f"Setting contacts with new current_values: {current_values}")
        for contact in self._all_contacts:
            if contact.name in current_values:
                contact.current = current_values[contact.name]

    @property
    def voltages(self) -> dict:
        """Returns the voltage values of each contact.

        Returns
        -------
        dict
        """
        return {contact.name: contact.voltage for contact in self._all_contacts}

    @voltages.setter
    def voltages(self, voltage_values: dict) -> None:
        """Set voltage value contacts.

        Parameters
        ----------
        voltage_values : dict
            Voltage values. Not all contacts have to be present in the
            dictionary.
        """
        _logger.debug(f"Setting contacts with new voltage_values: {voltage_values}")
        for contact in self._all_contacts:
            if contact.name in voltage_values:
                contact.voltage = voltage_values[contact.name]

    @property
    def surface_impedances(self) -> dict:
        """Returns the floating impedance values of each contact.

        Returns
        -------
        dict
        """
        return {
            contact.name: contact.surface_impedance for contact in self._all_contacts
        }

    @surface_impedances.setter
    def surface_impedances(self, impedance_values: dict) -> None:
        """Set surface impedance value contacts.

        Parameters
        ----------
        impedance_values : dict
            Surface impedance values. Not all contacts have to be present in the
            dictionary.
        """
        for contact in self._all_contacts:
            if contact.name in impedance_values:
                contact.impedance = impedance_values[contact.name]

    def __getitem__(self, name):
        for contact in self._all_contacts:
            if name == contact.name:
                return contact

    def __iter__(self):
        return iter(self._all_contacts)
