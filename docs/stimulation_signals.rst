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

Supported signal classes
------------------------

The package currently includes rectangular, trapezoidal, and triangular signal
classes. In practice, rectangular DBS pulses are the most mature and are the
best starting point for external users.

Typical settings
----------------

Common signal-related settings include:

- signal type
- stimulation frequency
- pulse width
- counter-pulse or inter-pulse settings when applicable
- whether the stimulation is current-controlled or voltage-controlled

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
