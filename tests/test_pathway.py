import h5py
import numpy as np
import pytest

from ossdbs.point_analysis.pathway import Pathway
from ossdbs.point_analysis.time_results import TimeResult


def create_one_pop_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        grp1 = f.create_group("pop1")
        ax0 = grp1.create_dataset("axon0", data=np.array([[0, 0, 0], [1, 1, 1]]))
        ax0.attrs["inx"] = 0
        ax1 = grp1.create_dataset("axon1", data=np.array([[2, 2, 2], [300, 3, 3]]))
        ax1.attrs["inx"] = 1


def create_two_pop_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        # First population: 2 axons
        grp1 = f.create_group("pop1")
        ax0_grp1 = grp1.create_dataset("axon0", data=np.array([[0, 0, 0], [1, 1, 1]]))
        ax0_grp1.attrs["inx"] = 0
        ax1_grp1 = grp1.create_dataset("axon1", data=np.array([[2, 2, 2], [3, 3, 3]]))
        ax1_grp1.attrs["inx"] = 1

        # Second population: 3 axons
        grp2 = f.create_group("pop2")
        ax0_grp2 = grp2.create_dataset("axon0", data=np.array([[4, 4, 4], [5, 5, 5]]))
        ax0_grp2.attrs["inx"] = 0
        ax1_grp2 = grp2.create_dataset("axon1", data=np.array([[6, 6, 6], [7, 7, 7]]))
        ax1_grp2.attrs["inx"] = 1
        ax2_grp2 = grp2.create_dataset("axon2", data=np.array([[8, 8, 8], [9, 9, 9]]))
        ax2_grp2.attrs["inx"] = 2


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


def create_three_axon_dummy_h5(filepath):
    with h5py.File(filepath, "w") as f:
        grp1 = f.create_group("pop1")
        for i in range(3):
            axon = grp1.create_dataset(
                f"axon{i}", data=np.array([[i, 0, 0], [i, 1, 1]], dtype=float)
            )
            axon.attrs["inx"] = i


@pytest.fixture
def h5_file_3_axons(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("h5data3")
    h5_file_3_axons = base_dir / "h5_file_3_axons.h5"
    create_three_axon_dummy_h5(h5_file_3_axons)
    return h5_file_3_axons


def export_potentials(pathway, potential, tmp_path, n_axons):
    """Run _write_file and read Potential[V] back per axon, None if not exported."""
    data = TimeResult(
        time_steps=np.array([0.0]),
        points=np.zeros((len(potential), 3)),
        inside_csf=np.zeros(len(potential), dtype=int),
        inside_encap=np.zeros(len(potential), dtype=int),
        potential=potential,
    )
    out_file = tmp_path / "exported.h5"
    with h5py.File(out_file, "w") as f:
        pathway._write_file(data, f)

    exported = []
    with h5py.File(out_file, "r") as f:
        for i in range(n_axons):
            group = f["pop1"][f"axon{i}"]
            if "Potential[V]" not in group:
                exported.append(None)
            else:
                exported.append(float(np.asarray(group["Potential[V]"])[0, 0]))
    return exported


def test_encap_marked_axon_does_not_shift_potential_export(h5_file_3_axons, tmp_path):
    """A clean axon must receive its own potential, not an earlier axon's.

    filter_csf_encap marks axons without removing their points from the lattice, so an
    encap-marked axon still owns its rows in data.potential. _create_datasets gates both
    the Potential[V] export and the advance of the `start` cursor on the same condition,
    `axon.status != -1`, so those rows are never consumed and every later axon is shifted.
    """
    pathway = Pathway(h5_file_3_axons)

    # One row per lattice point; axon k owns rows 2k and 2k+1 and should export 10*k.
    potential = np.repeat(np.arange(3) * 10.0, 2).reshape(6, 1)

    inside_csf = np.zeros(6, dtype=int)
    inside_encap = np.zeros(6, dtype=int)
    inside_encap[0] = 1  # axon0 touches the encapsulation layer
    pathway.filter_csf_encap(inside_csf, inside_encap)

    exported = export_potentials(pathway, potential, tmp_path, 3)

    assert exported[0] is None  # marked, deliberately not exported
    assert exported[1] == 10.0
    assert exported[2] == 20.0


def test_geometry_and_encap_marks_are_not_interchangeable(h5_file_3_axons, tmp_path):
    """The two sources of status == -1 need opposite cursor handling.

    axon0 lies outside the domain, so filter_for_geometry drops its points and it owns no
    rows: the cursor must NOT advance for it. axon1 lies inside the encapsulation layer,
    so its rows are retained and the cursor MUST advance. Advancing unconditionally would
    fix the second case and break the first.
    """
    pathway = Pathway(h5_file_3_axons)

    points = np.array([[i, j, j] for i in range(3) for j in range(2)], dtype=float)
    mask = np.zeros_like(points, dtype=bool)
    mask[0:2, :] = True  # axon0 lies outside the domain
    lattice = pathway.filter_for_geometry(np.ma.MaskedArray(points, mask=mask))
    assert lattice.shape == (4, 3)

    inside_csf = np.zeros(4, dtype=int)
    inside_encap = np.zeros(4, dtype=int)
    inside_encap[0] = 1  # axon1's first retained point
    pathway.filter_csf_encap(inside_csf, inside_encap)

    # The lattice holds axon1's two rows, then axon2's two rows.
    potential = np.repeat([10.0, 20.0], 2).reshape(4, 1)
    exported = export_potentials(pathway, potential, tmp_path, 3)

    assert exported[0] is None  # outside the domain: owns no rows
    assert exported[1] is None  # encap-marked, but its rows ARE in the lattice
    assert exported[2] == 20.0


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
