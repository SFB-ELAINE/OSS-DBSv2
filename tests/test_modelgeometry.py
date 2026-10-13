import logging

import numpy as np
import pytest

import ossdbs
from ossdbs.utils.settings import Settings


class TestModelGeometry:
    """Class for testing model geometry."""

    @pytest.fixture
    def modelGeometry(self):
        input_settings = {
            "BrainRegion": {
                "Center": {"x[mm]": 5.0, "y[mm]": 0.0, "z[mm]": 0.0},
                "Dimension": {"x[mm]": 50.0, "y[mm]": 50.0, "z[mm]": 50.0},
                "Shape": "Box",
            },
            "Electrodes": [
                {
                    "Name": "AbbottStJudeActiveTip6142_6145",
                    "Rotation[Degrees]": 0.0,
                    "Direction": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 1.0},
                    "TipPosition": {"x[mm]": 0.0, "y[mm]": 0.0, "z[mm]": 0.0},
                    "Contacts": [
                        {
                            "Contact_ID": 1,
                            "Active": True,
                            "Current[A]": 0.0,
                            "Voltage[V]": 1.0,
                            "Floating": False,
                            "SurfaceImpedance": {"Model": None, "Parameters": {}},
                            "MaxMeshSizeEdge": 0.01,
                        },
                        {
                            "Contact_ID": 2,
                            "Active": True,
                            "Current[A]": 0.0,
                            "Voltage[V]": 0.0,
                            "Floating": False,
                            "SurfaceImpedance": {"Model": None, "Parameters": {}},
                            "MaxMeshSizeEdge": 0.01,
                        },
                    ],
                    "EncapsulationLayer": {
                        "Thickness[mm]": 0.1,
                        "Material": "Blood",
                        "DielectricModel": "ColeCole4",
                        "MaxMeshSize": 0.5,
                    },
                },
                {
                    "Name": "AbbottStJudeActiveTip6142_6145",
                    "Rotation[Degrees]": 0,
                    "Direction": {"x[mm]": 0, "y[mm]": 0, "z[mm]": 1},
                    "TipPosition": {"x[mm]": 10.0, "y[mm]": 0.0, "z[mm]": 0.0},
                    "EncapsulationLayer": {"Thickness[mm]": 0.1},
                },
            ],
            "Mesh": {
                "LoadMesh": False,
                "MeshElementOrder": 2,
                "MeshingHypothesis": {"Type": "Default", "MaxMeshSize": 10.0},
                "MeshSize": {
                    "Edges": {},
                    "Faces": {"E1C1": 0.1},
                    "Volumes": {"Brain": 0.5},
                },
            },
        }

        settings = Settings(input_settings).complete_settings()

        electrodes = ossdbs.generate_electrodes(settings)
        brain_region = ossdbs.create_bounding_box(settings["BrainRegion"])
        shape = settings["BrainRegion"]["Shape"]
        brain = ossdbs.BrainGeometry(shape, brain_region)
        geometry = ossdbs.ModelGeometry(brain, electrodes)
        ossdbs.set_contact_and_encapsulation_layer_properties(settings, geometry)

        return geometry, settings, brain_region, electrodes

    def test_geometry(self, modelGeometry):
        """Test geometry()."""
        settings = modelGeometry[1]
        brain_region = modelGeometry[2]
        dimension = settings["BrainRegion"]["Dimension"]
        brain_shape = settings["BrainRegion"]["Shape"]
        shape = (dimension["x[mm]"], dimension["y[mm]"], dimension["z[mm]"])
        x, y, z = np.subtract(brain_region.end, brain_region.start) / 2

        brain_vol = None
        if brain_shape == "Box":
            brain_vol = shape[0] * shape[1] * shape[2]
        elif brain_shape == "Sphere":
            radius = np.min([x, y, z])
            brain_vol = 4 / 3 * np.pi * radius**3
        elif brain_shape == "Ellipsoid":
            brain_vol = 4 / 3 * np.pi * x * y * z

        electrode_vol = 0
        electrodes = modelGeometry[3]
        for electrode in electrodes:
            lead_radius = electrode._parameters.lead_diameter * 0.5
            total_length = np.max([x, y, z])
            height = total_length - lead_radius
            electrode_vol += (np.pi * lead_radius**2 * height) + (
                4 / 3 * np.pi * lead_radius**3 * 0.5
            )

        desired = brain_vol - electrode_vol
        actual = modelGeometry[0].geometry.shape.mass
        tolerance = 1e-5
        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def test_get_contact_index(self, modelGeometry):
        """Test get_contact_index()."""
        geometry = modelGeometry[0]
        contact_name = {contact.name for contact in geometry.contacts}

        actual = {geometry.get_contact_index(contact) for contact in contact_name}
        desired = set(range(len(geometry.contacts)))

        np.testing.assert_equal(actual, desired)

    def test_update_contact(self, modelGeometry):
        """Test update_contact()."""
        new_properties = {
            "Active": True,
            "Current[A]": 2.0,
            "Floating": False,
            "Voltage[V]": 4.0,
            "SurfaceImpedance": {"Model": "CPE_dl", "Parameters": {"dl_k": 1.5e6}},
        }
        geometry = modelGeometry[0]
        geometry.update_contact(0, new_properties)

        contact = geometry.contacts[0]
        assert contact.active is True
        assert contact.current == 2.0
        assert contact.floating is False
        assert contact.voltage == 4.0
        assert contact.surface_impedance_model == "CPE_dl"
        assert contact.surface_impedance_parameters == {"dl_k": 1.5e6}

    def test_get_encapsulation_layer_index(self, modelGeometry):
        """Test encapsulation_layer_index()."""
        geometry = modelGeometry[0]
        layer_name = {layer.name for layer in geometry.encapsulation_layers}

        actual = {geometry.get_encapsulation_layer_index(layer) for layer in layer_name}
        desired = set(range(len(geometry.encapsulation_layers)))

        np.testing.assert_equal(actual, desired)

    def test_update_encapsulation_layer(self, modelGeometry):
        """Test update_encapsulation_layer()."""
        new_properties = {
            "Material": "Gray matter",
            "DielectricModel": "ColeCole3",
            "MaxMeshSize": 0.8,
        }
        geometry = modelGeometry[0]
        geometry.update_encapsulation_layer(0, new_properties)

        desired = set(new_properties.values())
        actual = set()
        actual.add(geometry.encapsulation_layers[0].material)
        actual.add(geometry.encapsulation_layers[0].dielectric_model)
        actual.add(geometry.encapsulation_layers[0].max_h)

        assert actual == desired

    def test_set_edge_mesh_sizes(self, modelGeometry):
        """Test set_edge_mesh_sizes()."""
        geometry = modelGeometry[0]
        test_val = 0.002
        # First edge whose name is not None
        test_edge = next(
            edge.name for edge in geometry._shape.edges if edge.name is not None
        )
        geometry.set_edge_mesh_sizes({test_edge: test_val})

        count = sum(1 for edge in geometry._shape.edges if edge.name == test_edge)
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [edge.maxh for edge in geometry._shape.edges if edge.name == test_edge]
        )

        return np.testing.assert_equal(actual, desired)

    def test_edge_mesh_sizes_affect_mesh(self, modelGeometry):
        """Edge mesh sizes must produce a finer mesh."""
        from ossdbs.fem import Mesh

        _, settings, brain_region, electrodes = modelGeometry
        mesh_params = settings["Mesh"]
        shape = settings["BrainRegion"]["Shape"]
        brain = ossdbs.BrainGeometry(shape, brain_region)

        # Reference mesh without edge refinement
        geo_ref = ossdbs.ModelGeometry(brain, electrodes)
        m_ref = Mesh(geo_ref.geometry, order=2)
        m_ref.generate_mesh(mesh_params)

        # Refined mesh: set small maxh on a contact edge
        brain2 = ossdbs.BrainGeometry(shape, brain_region)
        geo_fine = ossdbs.ModelGeometry(brain2, electrodes)
        test_edge = next(e.name for e in geo_fine._shape.edges if e.name is not None)
        geo_fine.set_edge_mesh_sizes({test_edge: 0.01})
        m_fine = Mesh(geo_fine.geometry, order=2)
        m_fine.generate_mesh(mesh_params)

        assert m_fine.ngsolvemesh.ne > m_ref.ngsolvemesh.ne

    def test_set_face_mesh_sizes(self, modelGeometry):
        """Test set_face_mesh_sizes()."""
        geometry = modelGeometry[0]
        test_val = 0.2
        test_face = next(
            face.name for face in geometry._shape.faces if face.name is not None
        )
        geometry.set_face_mesh_sizes({test_face: test_val})

        count = sum(1 for face in geometry._shape.faces if face.name == test_face)
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [face.maxh for face in geometry._shape.faces if face.name == test_face]
        )

        return np.testing.assert_equal(actual, desired)

    class _FakeFace:
        def __init__(self, name):
            self.name = name

    class _FakeBrainGeo:
        def __init__(self, face_names):
            self.faces = [TestModelGeometry._FakeFace(name) for name in face_names]

    class _FakeElectrode:
        def __init__(
            self,
            index,
            n_contacts,
            require_all_contacts=True,
            required_contact_indices=frozenset(),
        ):
            self.index = index
            self.n_contacts = n_contacts
            self.require_all_contacts = require_all_contacts
            self.required_contact_indices = required_contact_indices

    @pytest.fixture
    def bare_geometry(self):
        """A ModelGeometry instance without running __init__.

        check_brain_geo() only uses get_contact_name() (pure), the
        electrode/brain_geo arguments, and
        self._contacts_missing_from_geometry, so no real geometry is needed.
        """
        geometry = object.__new__(ossdbs.ModelGeometry)
        geometry._contacts_missing_from_geometry = set()
        return geometry

    def test_check_brain_geo_all_contacts_present(self, bare_geometry):
        """All contacts found: passes regardless of require_all_contacts."""
        electrode = self._FakeElectrode(index=1, n_contacts=2)
        brain_geo = self._FakeBrainGeo(["E1C1", "E1C2"])
        assert bare_geometry.check_brain_geo(brain_geo, electrode) is True
        assert bare_geometry._contacts_missing_from_geometry == set()

    def test_check_brain_geo_missing_contact_strict_fails(self, bare_geometry):
        """Default (require_all_contacts=True): any missing contact fails."""
        electrode = self._FakeElectrode(
            index=1, n_contacts=2, require_all_contacts=True
        )
        brain_geo = self._FakeBrainGeo(["E1C1"])
        assert bare_geometry.check_brain_geo(brain_geo, electrode) is False

    def test_check_brain_geo_missing_contact_relaxed_warns(self, bare_geometry, caplog):
        """require_all_contacts=False: a missing contact only warns."""
        electrode = self._FakeElectrode(
            index=1, n_contacts=2, require_all_contacts=False
        )
        brain_geo = self._FakeBrainGeo(["E1C1"])
        with caplog.at_level(logging.WARNING):
            assert bare_geometry.check_brain_geo(brain_geo, electrode) is True
        assert "E1C2" in caplog.text
        assert bare_geometry._contacts_missing_from_geometry == {"E1C2"}

    def test_check_brain_geo_all_missing_relaxed_still_fails(self, bare_geometry):
        """require_all_contacts=False: zero contacts present still fails."""
        electrode = self._FakeElectrode(
            index=1, n_contacts=2, require_all_contacts=False
        )
        brain_geo = self._FakeBrainGeo([])
        assert bare_geometry.check_brain_geo(brain_geo, electrode) is False

    def test_check_brain_geo_missing_required_contact_still_fails(self, bare_geometry):
        """require_all_contacts=False: a missing Active/Floating contact
        still fails, since silently dropping it would produce wrong FEM
        results rather than a clear error.
        """
        electrode = self._FakeElectrode(
            index=1,
            n_contacts=2,
            require_all_contacts=False,
            required_contact_indices={2},
        )
        brain_geo = self._FakeBrainGeo(["E1C1"])
        assert bare_geometry.check_brain_geo(brain_geo, electrode) is False

    def test_check_brain_geo_missing_non_required_contact_warns(
        self, bare_geometry, caplog
    ):
        """require_all_contacts=False: a missing contact that is neither
        Active nor Floating only warns, even when other contacts on the
        same electrode are required.
        """
        electrode = self._FakeElectrode(
            index=1,
            n_contacts=2,
            require_all_contacts=False,
            required_contact_indices={1},
        )
        brain_geo = self._FakeBrainGeo(["E1C1"])
        with caplog.at_level(logging.WARNING):
            assert bare_geometry.check_brain_geo(brain_geo, electrode) is True
        assert "E1C2" in caplog.text

    def test_update_contact_areas_missing_contact_strict_raises(self, bare_geometry):
        """A contact absent from the shape and not marked as allowed-missing
        still raises, matching the pre-existing strict behaviour.
        """
        bare_geometry._shape = self._FakeBrainGeo(["E1C1"])
        bare_geometry._contacts = [ossdbs.model_geometry.Contact(name="E1C2")]
        with pytest.raises(RuntimeError, match="Area for E1C2 not set"):
            bare_geometry.update_contact_areas()

    def test_update_contact_areas_missing_contact_relaxed_warns(
        self, bare_geometry, caplog
    ):
        """A contact recorded by check_brain_geo as allowed-missing only
        warns instead of raising.
        """
        bare_geometry._shape = self._FakeBrainGeo(["E1C1"])
        bare_geometry._contacts = [ossdbs.model_geometry.Contact(name="E1C2")]
        bare_geometry._contacts_missing_from_geometry = {"E1C2"}
        with caplog.at_level(logging.WARNING):
            bare_geometry.update_contact_areas()
        assert "E1C2" in caplog.text
        assert bare_geometry._contacts[0].area is None

    def test_set_volume_mesh_sizes(self, modelGeometry):
        """Test set_volume_mesh_sizes()."""
        geometry = modelGeometry[0]
        test_val = 1.2
        test_volume = next(
            solid.name for solid in geometry._shape.solids if solid.name is not None
        )
        geometry.set_volume_mesh_sizes({test_volume: test_val})

        count = sum(1 for solid in geometry._shape.solids if solid.name == test_volume)
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [
                solid.maxh
                for solid in geometry._shape.solids
                if solid.name == test_volume
            ]
        )

        return np.testing.assert_equal(actual, desired)
