

from ossdbs.spectrum_modes.spectrum_point_analysis import SpectrumMode
from ossdbs.spectrum_modes.spectrum_point_analysis import FullSpectrum, OctaveBandMode
from ossdbs.spectrum_modes.voltage_setting import VoltageSetterCurrentControlled
from ossdbs.spectrum_modes.voltage_setting import VoltageSetterVoltageControlled
from ossdbs.spectrum_modes.spectrum_point_analysis_bipolar_current_control \
    import FullSpectrum as FullSpectrumBiPolar
from ossdbs.spectrum_modes.spectrum_point_analysis_bipolar_current_control \
    import OctaveBandMode as OctaveBandBiPolar


class SpectrumFactory:

    def create(spectrum_mode: str,
               current_controlled: bool,
               n_active_contacts: int) -> SpectrumMode:

        voltage_setter = VoltageSetterVoltageControlled

        if current_controlled and n_active_contacts > 2:
            voltage_setter = VoltageSetterCurrentControlled

        if current_controlled and n_active_contacts < 3:
            if spectrum_mode == 'OctaveBand':
                return OctaveBandBiPolar(voltage_setter)
            return FullSpectrumBiPolar(voltage_setter)

        if spectrum_mode == 'OctaveBand':
            return OctaveBandMode(voltage_setter)
        return FullSpectrum(voltage_setter)
