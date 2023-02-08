import numpy as np
from dataclasses import dataclass


@dataclass
class OctaveBand:

    SQRT2 = np.sqrt(2)

    def __init__(self, frequency: float) -> None:
        self.frequency = frequency

    def lower_limit(self):
        return self.frequency / self.SQRT2

    def upper_limit(self):
        return self.frequency * self.SQRT2


frequencies = np.arange(5001) * 130
n_octaves = int(np.log2(len(frequencies) - 1))
octave_indices = 2 ** np.arange(0, n_octaves)
octave_frequencies = frequencies[octave_indices]
octave_bands = [OctaveBand(frequency) for frequency in octave_frequencies]
last_freq = octave_frequencies[-1] * 2
extended_bands = [OctaveBand(0)] + octave_bands + [OctaveBand(last_freq)]

data = np.zeros(len(frequencies))
for octave in extended_bands:
    lower_idx = int(octave.lower_limit() / 130 + 1)
    upper_idx = int(octave.upper_limit() / 130 + 1)
    data[lower_idx:upper_idx] = octave.center_freq

    print(len(frequencies[lower_idx:upper_idx]))

    print(frequencies[lower_idx], octave.center_freq)
print(frequencies[-1])

print(data.shape)