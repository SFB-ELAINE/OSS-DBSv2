
from typing import List
from ossdbs.electrodes.contacts import Contacts
from ossdbs.impedance_analysis.impedances import Impedances
from ossdbs.impedance_analysis.spectrum_modes.spectrum_mode import SpectrumMode
from ossdbs.stimmulation_signals.trapzoid_signal import Signal
from ossdbs.fem import VolumeConductor
import ngsolve
import numpy as np


class OctaveBand(SpectrumMode):

    class OctaveBand:

        SQRT2 = np.sqrt(2)

        def __init__(self, frequency: float) -> None:
            self.frequency = frequency

        def lower_limit(self) -> float:
            return self.frequency / self.SQRT2

        def upper_limit(self) -> float:
            return self.frequency * self.SQRT2

    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                contacts: Contacts
                ) -> Impedances:

        frequencies = signal.fft_frequncies()
        settings = self.__contact_settings(contacts)
        n_octaves = int(np.log2(len(frequencies) - 1)) + 1
        impedances = np.zeros((n_octaves + 1, len(settings)), dtype=complex)

        for index, contact_setting in enumerate(settings):
            impedances[0, index] = self.__compute_impedance(volume_conductor,
                                                            0.0,
                                                            contact_setting)

        octave_bands = self.__octave_bands(signal.frequency, n_octaves)

        for freq_index, octave_band in enumerate(octave_bands, 1):
            for index, contact_setting in enumerate(settings):
                impedance = self.__compute_impedance(volume_conductor,
                                                     octave_band.frequency,
                                                     contact_setting)
                impedances[freq_index, index] = impedance

        contact_sets = [[contact.name for contact in contacts_setting.active()]
                        for contacts_setting in settings]

        frequencies = [0] + list(signal.frequency * 2 ** np.arange(n_octaves))
        return Impedances(frequencies=frequencies,
                          imdedances=impedances,
                          contact_sets=contact_sets)

    def __octave_bands(self,
                       frequency: float,
                       n_octaves: int
                       ) -> List[OctaveBand]:
        octave_indices = 2 ** np.arange(0, n_octaves)
        octave_frequencies = frequency * octave_indices
        return [self.OctaveBand(freq) for freq in octave_frequencies]

    @staticmethod
    def __compute_impedance(volume_conductor: VolumeConductor,
                            frequency: float,
                            contacts: Contacts
                            ) -> float:
        solution = volume_conductor.compute_solution(frequency, contacts)
        field = ngsolve.grad(solution.potential)
        curr_dens_conj = ngsolve.Conj(solution.current_density)
        mesh = volume_conductor.mesh.ngsolvemesh()
        power = ngsolve.Integrate(field * curr_dens_conj, mesh)
        voltage = 1
        return voltage / power if power else 0

    def __contact_settings(self, contacts: Contacts):
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
