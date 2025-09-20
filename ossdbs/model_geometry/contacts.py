# Copyright 2023, 2024 Johannes Reding, Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

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
        * `area`: Will be set during the geometry creation process.

    Mesh related properties:

        * `max_h`: Maximum element size on surface
        * `edge_max_h`: Maximum element size on contact edges

    Volume conductor model related properties:

        * `active`: Whether it needs a Dirichlet BC.
        * `floating`: Whether the voltage shall be fixed but unknown.
        * `current`: Assigned or computed current value.
        * `voltage`: Assigned or computed voltage value.
        * `surface_impedance_model`: Assigned surface impedance model.
        * `surface_impedance_parameters`: Parameters for surface impedance model.
    """

    name: str
    area: Optional[float] = None
    max_h: float = 1e10  # Netgen default
    edge_max_h: float = 1e10
    active: bool = False
    floating: bool = False
    current: float = 0.0
    voltage: float = 0.0
    surface_impedance_model: Optional[str] = None
    surface_impedance_parameters: Optional[dict] = None

    def __str__(self):
        """Write properties to string."""
        contact_str = self.name
        contact_str += f"\nMax h: {self.max_h}"
        contact_str += f"\nEdge max h: {self.edge_max_h}"
        contact_str += f"\nActive: {self.active}"
        contact_str += f"\nFloating: {self.floating}"
        contact_str += f"\nCurrent: {self.current}"
        contact_str += f"\nVoltage: {self.voltage}"
        if self.surface_impedance_model is not None:
            contact_str += (
                f"\nSurface impedance model: {self.surface_impedance_model}\n"
            )
            contact_str += (
                f"\nSurface impedance parameters: {self.surface_impedance_parameters}\n"
            )
        return contact_str

    def get_surface_impedance(
        self, frequency: float, is_complex: bool
    ) -> float | complex:
        """Return surface impedance at fixed frequency."""
        try:
            import impedancefitter as ifit
        except ImportError as err:
            raise ImportError(
                "Please install impedancefitter to compute the surface impedance."
                "Ensure that impedancefitter and its dependencies were correctly "
                "installed."
            ) from err

        # TODO cache this variable
        ecm = ifit.get_equivalent_circuit_model(self.surface_impedance_model)
        # TODO add surface area to it?
        Z = ecm.eval(omega=2.0 * np.pi * frequency, **self.surface_impedance_parameters)
        surface_Z = complex(Z) * self.area
        if not is_complex:
            return np.abs(surface_Z)
        return surface_Z


def check_contact(contact: Contact):
    """Check if contact has a clear role."""
    if contact.active and contact.floating:
        raise ValueError(
            f"""The contact {contact.name} has multiple roles.
            Please make sure that contacts are either active,
            floating or none of the two."""
        )
    if contact.surface_impedance_model is not None:
        if contact.surface_impedance_parameters is None:
            raise ValueError(
                "Surface impedance model was provided without parameter dictionary."
            )
    else:
        if contact.surface_impedance_parameters is not None:
            _logger.warning(
                "Surface impedance parameter dictionary was provided without model. "
                "It will not be taken into account."
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
        self._update_contacts()

    def _update_contacts(self) -> None:
        """Update contacts after change."""
        # Dirichlet boundary conditions
        self._active = [contact for contact in self._all_contacts if contact.active]
        # Floating boundary conditions
        self._floating = [contact for contact in self._all_contacts if contact.floating]
        # Do not stimulate / Neumann boundary condition
        self._unused = [
            contact
            for contact in self._all_contacts
            if not contact.floating and not contact.active
        ]
        # if a contact has a surface impedance model,
        # set the surface impedance to active
        self._surface_impedance_active = False
        for contact in self._all_contacts:
            if contact.surface_impedance_model is not None:
                self._surface_impedance_active = True
                break

    def append(self, contact: Contact) -> None:
        """Add another contact."""
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
    def surface_impedance_active(self) -> bool:
        """Get status of surface impedance."""
        return self._surface_impedance_active

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

    def get_surface_impedances(self, frequency: float, is_complex: bool) -> dict:
        """Returns the floating impedance values of each contact at a fixed frequency.

        Returns
        -------
        dict
        """
        return {
            contact.name: contact.get_surface_impedance(frequency, is_complex)
            if contact.surface_impedance_model is not None
            else None
            for contact in self._all_contacts
        }

    def update_contact(self, name, floating=None, active=None):
        """Change type of contact."""
        contact = self.__getitem__(name)
        if floating is not None:
            contact.floating = floating
        if active is not None:
            contact.active = active
        self._update_contacts()

    def __getitem__(self, name):
        """Get contact by name."""
        for contact in self._all_contacts:
            if name == contact.name:
                return contact

    def __iter__(self):
        """Iterate over contacts."""
        return iter(self._all_contacts)

    def __str__(self):
        """Write info of all contacts to string."""
        contacts_str = ""
        for contact in self._all_contacts:
            contacts_str += str(contact)
        return contacts_str
