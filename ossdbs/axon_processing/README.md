# Axon models for pathway activation modelling

## Overview
NEURON (https://neuron.yale.edu/) is used to determine if an axon
is activated by the electric field that OSS-DBS computes.
Currently, two models are implemented: MRG2002 and McNeal1976.
Their NEURON files are in `neuron_templates/` and will be copied and
recompiled with every stimulation run.

### MRG2002
The MRG2002 comes with an involved axon model, which is implemented in 
`axon_default_MRG2002.py`.
The model parameters are only available for a limited number of fiber diameters. 

### McNeal1976
This model has a simpler structure and less parameters.

## Axon allocation

The axon morphology is implemented in `axon_models.py`.
There, the conversion from fibers is done.
Also, all toolsets to update the parameters when estimating
PAM or VTA are implemented.

## Neural stimulation

The NEURON interface is implemented in `neuron_model.py`.  
There, all axons can be tested w.r.t. their activation.
