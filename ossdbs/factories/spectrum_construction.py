

from ossdbs.spectrum_point_analysis import SpectrumMode
from ossdbs.spectrum_point_analysis import FullSpectrum, OctaveBandMode
from ossdbs.voltage_setting import VoltageSetterCurrentControlled
from ossdbs.voltage_setting import VoltageSetterVoltageControlled


class SpectrumFactory:

    def create(spectrum_mode: str, current_controlled: bool) -> SpectrumMode:

        voltage_setter = VoltageSetterVoltageControlled
        if current_controlled:
            voltage_setter = VoltageSetterCurrentControlled

        if spectrum_mode == 'OctaveBand':
            return OctaveBandMode(voltage_setter)
        return FullSpectrum(voltage_setter)
