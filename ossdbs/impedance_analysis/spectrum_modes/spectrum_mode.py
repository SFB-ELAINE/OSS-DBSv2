

from abc import ABC, abstractmethod
from ossdbs.electrodes.contacts import Contacts
from ossdbs.impedance_analysis.impedances import Impedances
from ossdbs.stimmulation_signals.trapzoid_signal import Signal
from ossdbs.fem import VolumeConductor


class SpectrumMode(ABC):

    @abstractmethod
    def compute(self,
                signal: Signal,
                volume_conductor: VolumeConductor,
                contacts: Contacts
                ) -> Impedances:
        pass
