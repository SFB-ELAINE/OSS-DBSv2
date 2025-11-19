"""Module to prepare and run pathway action models."""

from .axon_models import AxonModels
from .neuron_model import MRG2002, McNeal1976, NeuronSimulator
from .utilities import compare_pathways


def get_neuron_model(model_type: str, pathways_dict: dict, pathway_solution_dir: str):
    """Turn a string into a NEURON model."""
    if "MRG2002" in model_type:
        neuron_model = MRG2002(pathways_dict, pathway_solution_dir)
    elif "McNeal1976" in model_type:
        neuron_model = McNeal1976(pathways_dict, pathway_solution_dir)
    else:
        raise NotImplementedError(f"Model {model_type} not yet implemented.")
    return neuron_model


__all__ = ("MRG2002", "AxonModels", "McNeal1976", "NeuronSimulator", "compare_pathways")
