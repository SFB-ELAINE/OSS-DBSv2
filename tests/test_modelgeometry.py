import json
import os

import numpy as np
import pytest

import ossdbs
from ossdbs.utils.settings import Settings


class TestModelGeometry:
    """Class for testing model geometry."""

    @pytest.fixture
    def modelGeometry(self):
        FILE_PATH = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "test_modelgeometry.json"
        )
        with open(FILE_PATH) as json_file:
            input_settings = json.load(json_file)
        settings = Settings(input_settings).complete_settings()

        electrodes = ossdbs.generate_electrodes(settings)
        brain_region = ossdbs.create_bounding_box(settings["BrainRegion"])
        shape = settings["BrainRegion"]["Shape"]
        brain = ossdbs.BrainGeometry(shape, brain_region)

        return ossdbs.ModelGeometry(brain, electrodes)

    def test_get_contact_index(self, modelGeometry):
        """Test get_contact_index()."""
        contact_name = {contact.name for contact in modelGeometry.contacts}

        actual = {modelGeometry.get_contact_index(contact) for contact in contact_name}
        desired = set(range(len(modelGeometry.contacts)))

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
        modelGeometry.update_contact(0, new_properties)

        desired = {True, 2.0, False, 4.0, (1 + 1j)}
        actual = set()
        actual.add(modelGeometry.contacts[0].active)
        actual.add(modelGeometry.contacts[0].current)
        actual.add(modelGeometry.contacts[0].floating)
        actual.add(modelGeometry.contacts[0].voltage)
        actual.add(modelGeometry.contacts[0].surface_impedance)

        assert actual == desired

    def test_get_encapsulation_layer_index(self, modelGeometry):
        """Test encapsulation_layer_index()."""
        layer_name = {layer.name for layer in modelGeometry.encapsulation_layers}

        actual = {
            modelGeometry.get_encapsulation_layer_index(layer) for layer in layer_name
        }
        desired = set(range(len(modelGeometry.encapsulation_layers)))

        np.testing.assert_equal(actual, desired)

    def test_update_encapsulation_layer(self, modelGeometry):
        """Test update_encapsulation_layer()."""
        new_properties = {
            "Material": "Gray matter",
            "DielectricModel": "ColeCole3",
            "MaxMeshSize": 0.8,
        }
        modelGeometry.update_encapsulation_layer(0, new_properties)

        desired = set(new_properties.values())
        actual = set()
        actual.add(modelGeometry.encapsulation_layers[0].material)
        actual.add(modelGeometry.encapsulation_layers[0].dielectric_model)
        actual.add(modelGeometry.encapsulation_layers[0].max_h)

        assert actual == desired

    def test_set_edge_mesh_sizes(self, modelGeometry):
        """Test set_edge_mesh_sizes()."""
        test_val = 0.002
        test_edge = next(
            edge.name
            for edge in modelGeometry._geometry.shape.edges
            if edge.name is not None
        )
        modelGeometry.set_edge_mesh_sizes({test_edge: test_val})

        count = sum(
            1 for edge in modelGeometry._geometry.shape.edges if edge.name == test_edge
        )
        desired = {test_val for i in range(count)}
        actual = set()
        actual.update(
            edge.maxh
            for edge in modelGeometry._geometry.shape.edges
            if edge.name == test_edge
        )

        return np.testing.assert_equal(actual, desired)

    def test_set_face_mesh_sizes(self, modelGeometry):
        """Test set_face_mesh_sizes()."""
        test_val = 0.2
        test_face = next(
            face.name
            for face in modelGeometry._geometry.shape.faces
            if face.name is not None
        )
        modelGeometry.set_face_mesh_sizes({test_face: test_val})

        count = sum(
            1 for face in modelGeometry._geometry.shape.faces if face.name == test_face
        )
        desired = {test_val for i in range(count)}
        actual = set()
        actual.update(
            face.maxh
            for face in modelGeometry._geometry.shape.faces
            if face.name == test_face
        )

        return np.testing.assert_equal(actual, desired)

    def test_set_volume_mesh_sizes(self, modelGeometry):
        """Test set_volume_mesh_sizes()."""
        test_val = 1.2
        test_volume = next(
            solid.name
            for solid in modelGeometry._geometry.shape.solids
            if solid.name is not None
        )
        modelGeometry.set_volume_mesh_sizes({test_volume: test_val})

        count = sum(
            1
            for solid in modelGeometry._geometry.shape.solids
            if solid.name == test_volume
        )
        desired = {test_val for i in range(count)}
        actual = set()
        actual.update(
            solid.maxh
            for solid in modelGeometry._geometry.shape.solids
            if solid.name == test_volume
        )

        return np.testing.assert_equal(actual, desired)
