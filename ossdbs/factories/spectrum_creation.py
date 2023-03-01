

from ossdbs.point_analysis.spectrum import FullSpectrum, OctaveBandMode


class SpectrumFactory:

    def create(mode: str):
        if mode == 'OctaveBand':
            return OctaveBandMode()
        return FullSpectrum()
