import h5py
import numpy as np
import pytest

from ossdbs.point_analysis.pathway import Pathway


def create_one_pop_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        grp1 = f.create_group("pop1")
        grp1.create_dataset("axon0", data=np.array([[0, 0, 0], [1, 1, 1]]))
        grp1.create_dataset("axon1", data=np.array([[2, 2, 2], [300, 3, 3]]))


def create_two_pop_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        # First population: 2 axons
        grp1 = f.create_group("pop1")
        grp1.create_dataset("axon0", data=np.array([[0, 0, 0], [1, 1, 1]]))
        grp1.create_dataset("axon1", data=np.array([[2, 2, 2], [3, 3, 3]]))
        # Second population: 3 axons
        grp2 = f.create_group("pop2")
        grp2.create_dataset("axon0", data=np.array([[4, 4, 4], [5, 5, 5]]))
        grp2.create_dataset("axon1", data=np.array([[6, 6, 6], [7, 7, 7]]))
        grp2.create_dataset("axon2", data=np.array([[8, 8, 8], [9, 9, 9]]))


@pytest.fixture
def h5_file_1_pathway(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("h5data1")
    h5_file_1_pathway = base_dir / "h5_file_1_pathway.h5"
    create_one_pop_dummy_h5(h5_file_1_pathway)
    return h5_file_1_pathway


@pytest.fixture
def h5_file_2_pathway(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("h5data2")
    h5_file_2_pathway = base_dir / "h5_file_2_pathway.h5"
    create_two_pop_dummy_h5(h5_file_2_pathway)
    return h5_file_2_pathway


def test_pathway_init_and_population_names(h5_file_1_pathway):
    pathway = Pathway(h5_file_1_pathway)
    assert pathway.get_population_names() == ["pop1"]
    axon_names = pathway.get_axon_names()
    assert axon_names == [["axon0", "axon1"]]
    axon_numbers = pathway.get_axon_numbers()
    assert axon_numbers == [2]


def test_create_index_and_axon_length(h5_file_1_pathway):
    pathway = Pathway(h5_file_1_pathway)
    lattice = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
    index = pathway.create_index(lattice)
    assert index.shape == (4, 1)
    assert np.all(index[:2] == 0)
    assert np.all(index[2:] == 1)
    assert pathway.get_axon_length() == 2


def test_filter_for_geometry(h5_file_1_pathway):
    pathway = Pathway(h5_file_1_pathway)
    # Simulate a masked array: first axon inside, second axon outside
    data = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
    mask = np.array([[False], [False], [True], [True]])
    grid_pts = np.ma.MaskedArray(data, mask=np.broadcast_to(mask, data.shape))
    filtered = pathway.filter_for_geometry(grid_pts)
    # Only first axon should remain
    assert filtered.shape == (2, 3)
    assert np.allclose(filtered, np.array([[0, 0, 0], [1, 1, 1]]))


def test_filter_csf_encap(h5_file_1_pathway):
    pathway = Pathway(h5_file_1_pathway)
    # 4 points, first two are axon0, last two axon1
    inside_csf = np.array([0, 1, 0, 0])
    inside_encap = np.array([0, 0, 1, 0])
    pathway.filter_csf_encap(inside_csf, inside_encap)
    pop = pathway._populations[0]
    assert pop.axons[0].status == -2  # axon0: one point in csf
    assert pop.axons[1].status == -1  # axon1: one point in encap


def test_save_as_nifti(h5_file_1_pathway):
    pathway = Pathway(h5_file_1_pathway)
    with pytest.raises(NotImplementedError):
        pathway.save_as_nifti(np.array([1, 2, 3]), "dummy.nii")


def test_two_populations(h5_file_2_pathway):
    pathway = Pathway(h5_file_2_pathway)
    # Check two populations with correct number of axons exist
    assert pathway.get_population_names() == ["pop1", "pop2"]
    axon_names = pathway.get_axon_names()
    assert axon_names == [["axon0", "axon1"], ["axon0", "axon1", "axon2"]]
    axon_numbers = pathway.get_axon_numbers()
    assert axon_numbers == [2, 3]
