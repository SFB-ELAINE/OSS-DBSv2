

from ossdbs.impedance_analysis.spectrum import LogarithmScanning
from ossdbs.impedance_analysis.spectrum import OctaveBandMode
from ossdbs.impedance_analysis.spectrum import SpectrumMode


class SpectrumImpedanceFactory:

    def create(mode: str) -> SpectrumMode:
        if mode == 'OctaveBand':
            return OctaveBandMode()
        return LogarithmScanning()
