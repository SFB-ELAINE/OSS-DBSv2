Stimulation signals
===================

OSS-DBSv2 separates the stimulation waveform from the geometric and material
model. This makes it possible to study how different pulse definitions affect
the frequency-domain solve and the reconstructed time-domain results.

Signal concepts
---------------

The code distinguishes between:

- time-domain signal definitions such as rectangular, trapezoidal, or triangular
  pulses
- a frequency-domain representation used by the volume conductor solver

This separation is useful because the solver can evaluate the field at selected
frequencies and then reconstruct the signal in time when needed.

Supported signal types
----------------------

Time-domain pulse signals
^^^^^^^^^^^^^^^^^^^^^^^^^

These signals define a periodic waveform that is Fourier-transformed for the
frequency-domain FEM solve and optionally reconstructed in the time domain
afterwards.

- **Rectangle** (``"Type": "Rectangle"``) — standard rectangular DBS pulse.
  The most common and best-tested signal type.
- **Trapezoid** (``"Type": "Trapezoid"``) — trapezoidal pulse with a finite
  rise time controlled by ``PulseTopWidth[us]``.
- **Triangle** (``"Type": "Triangle"``) — triangular pulse (trapezoid with
  zero top width).

Common parameters for all time-domain signals:

- ``Frequency[Hz]`` — stimulation repetition rate (e.g. 130 Hz for standard
  DBS)
- ``PulseWidth[us]`` — duration of the primary pulse phase
- ``CounterPulseWidth[us]`` — duration of the charge-balancing counter pulse
  (0 to omit)
- ``InterPulseWidth[us]`` — gap between primary and counter pulse
- ``CounterAmplitude`` — amplitude of the counter pulse relative to the
  primary pulse (default 1.0)

Multisine
^^^^^^^^^

The multisine mode (``"Type": "Multisine"``) bypasses waveform generation
entirely. Instead, a list of discrete frequencies is solved with unit
amplitude at each:

.. code-block:: json

   "StimulationSignal": {
     "Type": "Multisine",
     "ListOfFrequencies": [130.0, 1000.0, 10000.0],
     "CurrentControlled": false
   }

This is useful for impedance spectroscopy, single-frequency studies, or when
the frequency content is known in advance. No time-domain reconstruction is
performed.

Frequency-domain settings
-------------------------

SpectrumMode
^^^^^^^^^^^^

For time-domain signals, ``SpectrumMode`` controls how many frequencies are
actually solved:

- ``"FullSpectrum"`` (default) — solves at every harmonic up to the cutoff
  frequency. Accurate but expensive for high cutoff frequencies.
- ``"OctaveBand"`` — solves only at octave-band centre frequencies and
  interpolates. Much faster with minimal loss of accuracy for typical DBS
  pulses.

CutoffFrequency
^^^^^^^^^^^^^^^^

``CutoffFrequency`` (default: 1e6 Hz) sets the upper limit of the Fourier
spectrum. Harmonics above this frequency are discarded. For most DBS
applications, 0.5–1 MHz is sufficient. Lower values reduce the number of
FEM solves.

CurrentControlled
^^^^^^^^^^^^^^^^^

``CurrentControlled`` selects between voltage-controlled and current-controlled
stimulation modes. See :ref:`stimulation-modes` in the volume conductor
documentation for a detailed description of all supported cases.

Related pages
-------------

- :doc:`volume_conductor_model` — stimulation modes, surface impedance, and
  solver configuration
- :doc:`input_settings` — ``StimulationSignal`` JSON settings
- :doc:`examples` — runnable example configurations

API reference
-------------

Signal classes
^^^^^^^^^^^^^^

.. automodule:: ossdbs.stimulation_signals.rectangle_signal
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.stimulation_signals.trapezoid_signal
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: ossdbs.stimulation_signals.triangle_signal
    :members:
    :undoc-members:
    :show-inheritance:

Base classes
^^^^^^^^^^^^

.. automodule:: ossdbs.stimulation_signals.signal
    :members:
    :undoc-members:
    :show-inheritance:

Helper functions
^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.stimulation_signals.utilities
    :members:
    :undoc-members:
