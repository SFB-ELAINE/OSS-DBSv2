

from ..spectrum_modes import SpectrumMode
from ..spectrum_modes import FullSpectrum, OctaveBand
from ..spectrum_modes import FullSpectrumBipolar
from ..spectrum_modes import OctaveBandBipolar
from ossdbs.point_analysis.voltage_setting import VoltageSetterCurrentControlled
from ossdbs.point_analysis.voltage_setting import VoltageSetterVoltageControlled


class SpectrumFactory:

    def create(spectrum_mode: str,
               current_controlled: bool,
               n_active_contacts: int) -> SpectrumMode:

        voltage_setter = VoltageSetterVoltageControlled

        if current_controlled and n_active_contacts > 2:
            voltage_setter = VoltageSetterCurrentControlled

        if current_controlled and n_active_contacts < 3:
            if spectrum_mode == 'OctaveBand':
                return OctaveBandBipolar(voltage_setter)
            return FullSpectrumBipolar(voltage_setter)

        if spectrum_mode == 'OctaveBand':
            return OctaveBand(voltage_setter)
        return FullSpectrum(voltage_setter)
