

from ossdbs.impedance_analysis.spectrum import LogarithmScanning, OctaveBandMode


class SpectrumImpedanceFactory:

    def create(mode: str):
        if mode == 'OctaveBand':
            return OctaveBandMode()
        return LogarithmScanning()
