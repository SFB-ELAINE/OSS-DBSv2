
from .signal import TimeDomainSignal
import numpy as np


class TrapezoidSignal(TimeDomainSignal):
    """Represents trapezoid signal.

    Parameters
    ----------
    frequency : float
        Frequency [Hz] of the signal.
    pulse_width : float
        Relative pulse width of one period.
    top_width: float
        Relative top width of trapzoid signal.
    counter_pulse_width: float
        Relative width of counter pulse of one period.
    inter_pulse_width: float
        Relative width between pulse and counter pulse of one period.
    """
    def __init__(self,
                 frequency: float,
                 pulse_width: float,
                 top_width: float,
                 counter_pulse_width: float,
                 inter_pulse_width: float) -> None:
        self._top_width = top_width
        super().__init__(frequency, pulse_width, inter_pulse_width, counter_pulse_width)

    def _generate_samples(self, sample_spacing: float) -> np.ndarray:
        """Generate samples which follow the signal form.

        Parameters
        ----------
        sample_spacing : float
            Timestep [s] between two samples.

        Returns
        -------
        np.ndarray
            Samples for one period.
        """

        spacing = min(abs(sample_spacing), 1) * self._frequency
        n_samples = int(1 / spacing) if 0 < spacing < 1 else 1

        pulse_length = int(self._pulse_width * n_samples)
        pulse_top_length = int(self._top_width * n_samples)
        pulse = self.__create_pulse(pulse_length, pulse_top_length)

        space_length = int(self._space_width * n_samples)
        counter_length = int(self._counter_pulse_width * n_samples)
        f = self._top_width / self._pulse_width if self._pulse_width else 0
        counter_top_length = int(f * self._counter_pulse_width * n_samples)
        counter_value = -pulse_length / counter_length if counter_length else 0
        counter_pulse = counter_value * self.__create_pulse(counter_length,
                                                            counter_top_length)
        counter_start = pulse_length + space_length
        counter_end = counter_start + counter_length

        signal = np.zeros(n_samples)
        signal[:pulse_length] = pulse
        n_counter_samples = n_samples - counter_start
        signal[counter_start:counter_end] = counter_pulse[:n_counter_samples]
        return signal

    @staticmethod
    def __create_pulse(pulse_length: int, pulse_top_length: int) -> np.ndarray:

        if not pulse_top_length:
            return np.array([0] * pulse_length)

        ramp_length = int((pulse_length - pulse_top_length) * 0.5)
        step_size = 1 / (ramp_length + 1)
        ramp = np.arange(step_size, 1, step_size)[:ramp_length]
        return np.concatenate((ramp, [1] * pulse_top_length, np.flip(ramp)))
