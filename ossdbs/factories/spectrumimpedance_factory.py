

from ossdbs.spectrum_modes.spectrum_impedance_analysis import LogarithmScanning
from ossdbs.spectrum_modes.spectrum_impedance_analysis import OctaveBandMode
from ossdbs.spectrum_modes.spectrum_impedance_analysis import SpectrumMode


class SpectrumImpedanceFactory:

    def create(mode: str) -> SpectrumMode:
        if mode == 'OctaveBand':
            return OctaveBandMode()
        return LogarithmScanning()
