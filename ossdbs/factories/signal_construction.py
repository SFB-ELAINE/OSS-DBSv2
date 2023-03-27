

from ossdbs.stimulation_signal import Signal
from ossdbs.stimulation_signal import RectangleSignal
from ossdbs.stimulation_signal import TrapzoidSignal
from ossdbs.stimulation_signal import TriangleSignal


class SignalFactory:

    @staticmethod
    def create(signal_parameters: dict) -> Signal:

        parameters = signal_parameters
        frequency = parameters['Frequency[Hz]']
        pulse_width = parameters['PulseWidth[µs]'] * 1e-6 * frequency
        counter_width = parameters['CounterPulseWidth[µs]'] * 1e-6 * frequency
        inter_width = parameters['InterPulseWidth[µs]'] * 1e-6 * frequency
        top_width = parameters['PulseTopWidth[µs]'] * 1e-6 * frequency

        if parameters['Type'] == 'Trapzoid':
            return TrapzoidSignal(frequency=frequency,
                                  pulse_width=pulse_width,
                                  top_width=top_width,
                                  counter_pulse_width=counter_width,
                                  inter_pulse_width=inter_width)

        if parameters['Type'] == 'Triangle':
            return TriangleSignal(frequency=frequency,
                                  pulse_width=pulse_width,
                                  counter_pulse_width=counter_width,
                                  inter_pulse_width=inter_width)

        return RectangleSignal(frequency=frequency,
                               pulse_width=pulse_width,
                               counter_pulse_width=counter_width,
                               inter_pulse_width=inter_width)
