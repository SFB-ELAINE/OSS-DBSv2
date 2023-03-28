

from ossdbs.electrodes.contacts import Contact, Contacts
from ossdbs.electrodes.electrodes import Electrodes


class ContactsFactory:

    def __init__(self, electrodes: Electrodes):
        self.__electrodes = electrodes

    def create(self, case_grounding: dict) -> Contacts:
        name = 'BrainSurface'
        active = case_grounding['Active']
        floating = False
        current = case_grounding['Current[A]']
        voltage = case_grounding['Voltage[V]']
        surface_impedance = 0.0 + 0.0j
        contact = Contact(name=name,
                          active=active,
                          floating=floating,
                          current=current,
                          voltage=voltage,
                          surface_impedance=surface_impedance)

        contacts = self.__electrodes.contacts()
        contacts.append(contact)
        return contacts
