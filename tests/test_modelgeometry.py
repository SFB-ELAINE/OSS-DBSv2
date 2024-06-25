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
                            "SurfaceImpedance[Ohmm]": {"real": 0.0, "imag": 0.0},
                            "MaxMeshSizeEdge": 0.01,
                        },
                        {
                            "Contact_ID": 2,
                            "Active": True,
                            "Current[A]": 0.0,
                            "Voltage[V]": 0.0,
                            "Floating": False,
                            "SurfaceImpedance[Ohmm]": {"real": 0.0, "imag": 0.0},
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
            "SurfaceImpedance[Ohmm]": {"real": 1, "imag": 1},
        }
        geometry = modelGeometry[0]
        geometry.update_contact(0, new_properties)

        desired = {True, 2.0, False, 4.0, (1 + 1j)}
        actual = set()
        actual.add(geometry.contacts[0].active)
        actual.add(geometry.contacts[0].current)
        actual.add(geometry.contacts[0].floating)
        actual.add(geometry.contacts[0].voltage)
        actual.add(geometry.contacts[0].surface_impedance)

        assert actual == desired

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
            edge.name
            for edge in geometry._geometry.shape.edges
            if edge.name is not None
        )
        geometry.set_edge_mesh_sizes({test_edge: test_val})

        count = sum(
            1 for edge in geometry._geometry.shape.edges if edge.name == test_edge
        )
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [
                edge.maxh
                for edge in geometry._geometry.shape.edges
                if edge.name == test_edge
            ]
        )

        return np.testing.assert_equal(actual, desired)

    def test_set_face_mesh_sizes(self, modelGeometry):
        """Test set_face_mesh_sizes()."""
        geometry = modelGeometry[0]
        test_val = 0.2
        test_face = next(
            face.name
            for face in geometry._geometry.shape.faces
            if face.name is not None
        )
        geometry.set_face_mesh_sizes({test_face: test_val})

        count = sum(
            1 for face in geometry._geometry.shape.faces if face.name == test_face
        )
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [
                face.maxh
                for face in geometry._geometry.shape.faces
                if face.name == test_face
            ]
        )

        return np.testing.assert_equal(actual, desired)

    def test_set_volume_mesh_sizes(self, modelGeometry):
        """Test set_volume_mesh_sizes()."""
        geometry = modelGeometry[0]
        test_val = 1.2
        test_volume = next(
            solid.name
            for solid in geometry._geometry.shape.solids
            if solid.name is not None
        )
        geometry.set_volume_mesh_sizes({test_volume: test_val})

        count = sum(
            1 for solid in geometry._geometry.shape.solids if solid.name == test_volume
        )
        desired = np.array([test_val for i in range(count)])
        actual = np.array(
            [
                solid.maxh
                for solid in geometry._geometry.shape.solids
                if solid.name == test_volume
            ]
        )

        return np.testing.assert_equal(actual, desired)
