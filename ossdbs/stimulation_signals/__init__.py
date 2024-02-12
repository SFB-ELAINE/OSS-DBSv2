# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Stimulation signal utility."""

from .rectangle_signal import RectangleSignal
from .signal import FrequencyDomainSignal, TimeDomainSignal
from .trapezoid_signal import TrapezoidSignal
from .triangle_signal import TriangleSignal
from .utilities import (
    generate_signal,
    get_indices_in_octave_band,
    get_maximum_octave_band_index,
    get_minimum_octave_band_index,
    get_octave_band_indices,
    get_timesteps,
    reconstruct_time_signals,
    retrieve_time_domain_signal_from_fft,
)

__all__ = [
    "TimeDomainSignal",
    "FrequencyDomainSignal",
    "RectangleSignal",
    "TrapezoidSignal",
    "TriangleSignal",
    "generate_signal",
    "retrieve_time_domain_signal_from_fft",
    "reconstruct_time_signals",
    "get_octave_band_indices",
    "get_indices_in_octave_band",
    "get_minimum_octave_band_index",
    "get_maximum_octave_band_index",
    "get_timesteps",
]
