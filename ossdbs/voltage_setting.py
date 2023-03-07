
from ossdbs.contacts import Contacts
from ossdbs.volume_conductor import VolumeConductor
from abc import ABC, abstractmethod
from typing import List
import numpy as np
import ngsolve


class VoltageSetter(ABC):

    @abstractmethod
    def set_voltages(cls,
                     frequency: float,
                     contacts: Contacts,
                     volume_conductor: VolumeConductor
                     ) -> Contacts:
        pass


class VoltageSetterVoltageControlled(VoltageSetter):

    @classmethod
    def set_voltages(cls,
                     frequency: float,
                     contacts: Contacts,
                     volume_conductor: VolumeConductor
                     ) -> Contacts:
        return contacts.copy()


class VoltageSetterCurrentControlled(VoltageSetter):

    @classmethod
    def set_voltages(cls,
                     frequency: float,
                     contacts: Contacts,
                     volume_conductor: VolumeConductor
                     ) -> Contacts:
        currents = [cls.__currents(frequency, contacts, volume_conductor)
                    for contacts in cls.__contact_settings()]
        admitances = np.array(currents).T
        impedance_matrix = np.linalg.inv(admitances)
        active_contacts = contacts.active_contacts()
        currents = [contact.current for contact in active_contacts]
        voltages = np.dot(impedance_matrix, currents)
        contacts_copy = contacts.copy()
        for contact, voltage in zip(active_contacts, voltages):
            contacts_copy.set_voltage(contact.name, voltage)
        return contacts_copy

    @classmethod
    def __currents(cls,
                   volume_conductor: VolumeConductor,
                   contacts: Contacts,
                   frequency: float
                   ) -> List[float]:
        solution = volume_conductor.compute_solution(frequency, contacts)
        density = ngsolve.BoundaryFromVolumeCF(solution.current_density)
        ng_mesh = volume_conductor.mesh.ngsolvemesh()
        return [cls.__contact_current(ng_mesh, density, contact)
                for contact in contacts.active_contacts()]

    @staticmethod
    def __contact_current(ng_mesh: ngsolve.Mesh,
                          density: ngsolve.CoefficientFunction,
                          boundary: str
                          ) -> ngsolve.CoefficientFunction:
        region = ng_mesh.Boundaries(boundary)
        current_v = ngsolve.Integrate(density, ng_mesh, definedon=region)
        return ngsolve.specialcf.normal(current_v, 3)

    @staticmethod
    def __contact_settings(contacts: Contacts) -> List[Contacts]:
        default_contacts = contacts.copy()
        for contact in default_contacts.active():
            default_contacts.set_voltage(contact.name, 0)
        return [default_contacts.copy().set_voltage(contact.name, 1.0)
                for contact in default_contacts.active()]
