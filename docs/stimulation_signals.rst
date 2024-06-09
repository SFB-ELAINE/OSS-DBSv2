Stimulation signals
===================

Signal types
------------

There are different types of signals.
Currently, only mono- and biphasic rectangular pulses are tested/usable.

API reference
-------------

Rectangular pulses
^^^^^^^^^^^^^^^^^^

.. automodule:: ossdbs.stimulation_signals.RectangleSignal
    :members:
    :undoc-members:
    :show-inheritance:

Base classes
^^^^^^^^^^^^

.. automodule:: ossdbs.stimulation_signals.TimeDomainSignal
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: ossdbs.stimulation_signals.FrequencyDomainSignal

Helper functions
^^^^^^^^^^^^^^^^

.. automethod:: ossdbs.stimulation_signals.retrieve_time_domain_signal_from_fft
.. automethod:: ossdbs.stimulation_signals.reconstruct_time_signals
.. automethod:: ossdbs.stimulation_signals.get_octave_band_indices
.. automethod:: ossdbs.stimulation_signals.get_indices_in_octave_band
.. automethod:: ossdbs.stimulation_signals.get_minimum_octave_band_index
.. automethod:: ossdbs.stimulation_signals.get_maximum_octave_band_index
.. automethod:: ossdbs.stimulation_signals.get_timesteps
