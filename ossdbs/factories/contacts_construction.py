

from ossdbs.electrodes.contacts import Contact, Contacts
from ossdbs.electrodes.electrodes import Electrodes


class ContactsFactory:

    def __init__(self, electrodes: Electrodes):
        self.__electrodes = electrodes

    def create(self, case_grounding: dict) -> Contacts:
        name = 'BrainSurface'
        active = case_grounding['Active']
        floating = case_grounding['Floating'] and not active
        current = case_grounding['Current[A]']
        voltage = case_grounding['Voltage[V]']
        real = case_grounding['SurfaceImpedance[Ωm]']['real']
        imag = case_grounding['SurfaceImpedance[Ωm]']['imag']
        surface_impedance = real + 1j * imag
        contact = Contact(name=name,
                          active=active,
                          floating=floating,
                          current=current,
                          voltage=voltage,
                          surface_impedance=surface_impedance)

        if not contact.active:
            return self.__electrodes.contacts()

        contacts = self.__electrodes.contacts()
        return Contacts(contacts.active() + contacts.floating() + [contact])
