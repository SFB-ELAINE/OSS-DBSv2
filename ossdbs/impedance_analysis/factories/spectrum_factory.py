

from ..spectrum_modes import LogarithmScanning
from ..spectrum_modes import OctaveBand
from ..spectrum_modes import SpectrumMode


class SpectrumFactory:

    def create(mode: str) -> SpectrumMode:
        if mode == 'OctaveBand':
            return OctaveBand()
        return LogarithmScanning()
