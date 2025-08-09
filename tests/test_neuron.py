import h5py
import numpy as np
import pytest

from ossdbs.axon_processing.neuron_model import McNeal1976, NeuronSimulator

neuron = pytest.importorskip("neuron")
pathway_dict = {}
pathway_dict_McNeal = {
    "n_Ranvier": [1],
    "axon_diams": [5.0],
    "Axon_Model_Type": "McNeal1976",
    "Name_prepared_neuron_array": "Allocated_axons.h5",
    "Neuron_model_array_prepared": True,
    "N_seeded_neurons": [1],
    "N_orig_neurons": [1],
    "connectome_name": "TestConnectome",
}


def create_dummy_h5(filepath):
    td_solution = h5py.File(filepath, "w")
    mock_signal = np.array([0, 1, 1, 1, 0, 0, 0, 0])
    td_solution.create_dataset("TimeSteps[s]", data=np.arange(len(mock_signal)) * 1e-6)
    grp1 = td_solution.create_group("pathway_1")
    grp1.create_dataset("Status", data=np.zeros(1))
    sub_group = grp1.create_group("axon0")
    sub_group.create_dataset("Potential[V]", data=mock_signal)
    sub_group.create_dataset("Points[mm]", data=np.array([[1, 1, 1]]))
    print(td_solution["TimeSteps[s]"])
    print(td_solution["pathway_1"]["axon0"])
    return td_solution


@pytest.fixture
def h5_file_1_pathway(tmp_path):
    h5_file_1_pathway = tmp_path / "h5_file_1_pathway.h5"
    td_solution = create_dummy_h5(h5_file_1_pathway)
    return td_solution


def test_neuron_model(tmp_path):
    with pytest.raises(TypeError):
        _ = NeuronSimulator(pathway_dict, tmp_path)


def test_initialise_mc_neal_model(tmp_path):
    _ = McNeal1976(pathway_dict_McNeal, tmp_path)


def test_process_mc_neal_model(tmp_path, h5_file_1_pathway):
    model = McNeal1976(pathway_dict_McNeal, tmp_path)
    model.process_pathways(h5_file_1_pathway)
    h5_file_1_pathway.close()
