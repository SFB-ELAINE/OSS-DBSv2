
from typing import List
from ossdbs.electrodes.contacts import Contacts
from ossdbs.impedance_analysis.impedances import Impedances
from ossdbs.impedance_analysis.spectrum_modes.spectrum_mode import SpectrumMode
from ossdbs.stimmulation_signals.trapzoid_signal import Signal
from ossdbs.fem import VolumeConductor
import ngsolve
import numpy as np


class LogarithmScanning(SpectrumMode):

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                contacts: Contacts
                ) -> Impedances:
        frequencies = self.__frequencies(signal)
        mesh = volume_conductor.mesh.ngsolvemesh
        settings = self.__contact_settings(contacts)
        impedances = np.zeros((len(frequencies), len(settings)), dtype=complex)

        for index, frequency in enumerate(frequencies):
            for set_index, contacts_setting in enumerate(settings):
                solution = volume_conductor.compute_solution(frequency,
                                                             contacts_setting)
                field = ngsolve.grad(solution.potential)
                curr_dens_conj = ngsolve.Conj(solution.current_density)
                power = ngsolve.Integrate(field * curr_dens_conj, mesh)
                impedances[index, set_index] = 1 / power if power else 0
                print(impedances[index])

        contact_sets = [[contact.name for contact in contacts_setting.active()]
                        for contacts_setting in settings]

        return Impedances(frequencies=frequencies,
                          imdedances=impedances,
                          contact_sets=contact_sets)

    def __frequencies(self, signal: Signal) -> List[float]:
        n_fft_frequencies = len(signal.fft_frequncies())
        n_intervalls = int(np.log10(n_fft_frequencies)) + 1
        frequencies = [0.0]
        for index in range(n_intervalls):
            interval__freq = signal.frequency * 10 ** index * np.arange(1, 10)
            frequencies.extend(interval__freq)
        return frequencies

    @staticmethod
    def __contact_settings(contacts: Contacts):
        settings = []
        for index, contact_1 in enumerate(contacts.active()[:-1]):
            for contact_2 in contacts.active()[index + 1:]:
                contact_1.voltage = 1
                contact_2.voltage = 0
                new_contacts = [contact_1, contact_2] + contacts.floating()
                settings.append(Contacts(new_contacts))

        return [Contacts([contact_1, contact_2] + contacts.floating())
                for index, contact_1 in enumerate(contacts.active()[:-1])
                for contact_2 in contacts.active()[index + 1:]]
