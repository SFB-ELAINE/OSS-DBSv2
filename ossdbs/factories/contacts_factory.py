

from ossdbs.electrodes.contacts import Contact


class CaseGroundContactFactory:

    @classmethod
    def create(cls, case_grounding: dict) -> Contact:
        name = 'BrainSurface'
        active = case_grounding['Active']
        floating = False
        current = case_grounding['Current[A]']
        voltage = case_grounding['Voltage[V]']
        surface_impedance = 0.0 + 0.0j
        return Contact(name=name,
                       active=active,
                       floating=floating,
                       current=current,
                       voltage=voltage,
                       surface_impedance=surface_impedance)
