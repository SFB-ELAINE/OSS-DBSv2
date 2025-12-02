import os
from typing import ClassVar

import ngsolve
import numpy as np

from ossdbs import generate_model_geometry

from .reference_directed_electrodes.tested_directions import (
    base_settings,
    directions,
    get_direction_dict,
)


def get_reference_and_model_geo(electrode_name, idx):
    direction = directions[idx]
    settings = base_settings.copy()
    settings["Electrodes"][0]["Name"] = electrode_name
    settings["Electrodes"][0]["Direction"] = get_direction_dict(direction)
    model_geo = generate_model_geometry(settings)
    mesh = ngsolve.Mesh(
        os.path.join(
            "tests",
            "electrode_tests",
            "reference_directed_electrodes",
            f"mesh_{electrode_name}_direction_{idx}.vol.gz",
        )
    )
    ref_geo = mesh.ngmesh.GetGeometry().shape
    return ref_geo, model_geo


# TODO is there a nicer way to not need to copy-paste
# the tests into each electrode file?
class TestElectrode:
    """Class for testing electrodes."""

    AbbottStJudeDirected: ClassVar[list[str]] = [
        "AbbottStJudeDirected6172",
        "AbbottStJudeDirected6173",
    ]
    Medtronic: ClassVar[list[str]] = ["Medtronic3387", "Medtronic3389", "Medtronic3391"]
    MedtronicSenSight: ClassVar[list[str]] = [
        "MedtronicSenSightB33005",
        "MedtronicSenSightB33015",
    ]
    NeuroPace: ClassVar[list[str]] = ["NeuroPaceDL344_3_5", "NeuroPaceDL344_10"]
    NeuroNexus: ClassVar[list[str]] = ["NeuroNexusA1x16_5mm_50_177"]
    PINSMedical: ClassVar[list[str]] = [
        "PINSMedicalL301",
        "PINSMedicalL302",
        "PINSMedicalL303",
    ]
    Dixi: ClassVar[list[str]] = [
        "DixiSEEG5",
        "DixiSEEG8",
        "DixiSEEG10",
        "DixiSEEG12",
        "DixiSEEG15",
        "DixiSEEG18",
    ]

    PMTsEEG2102: ClassVar[list[str]] = [
        "PMTsEEG2102_08",
        "PMTsEEG2102_10",
        "PMTsEEG2102_12",
        "PMTsEEG2102_14",
        "PMTsEEG2102_16",
    ]
    SceneRay1242: ClassVar[list[str]] = ["SceneRay1242"]

    def check_rename_boundaries(self, electrode, electrode_name):
        """Check whether set_contact_names() works."""
        changes = {
            "Body": "RenamedBody",
            "Contact_1": "RenamedContact_1",
            "NonExistingPart": "NonExistingPart",
        }

        electrode.set_contact_names(changes)
        n_contacts = electrode.n_contacts
        geometry = electrode.geometry
        faces = [face.name for face in geometry.faces]

        desired = {"RenamedBody", "RenamedContact_1"}

        desired.update({f"Contact_{i + 2}" for i in range(n_contacts - 1)})

        assert desired == set(faces)

    def check_contacts(self, electrode, electrode_name):
        """Check the number and names of contacts."""
        n_contacts = electrode.n_contacts
        geometry = electrode.geometry
        faces = [face.name for face in geometry.faces]

        desired = {"Body"}
        desired.update({f"Contact_{i + 1}" for i in range(n_contacts)})

        # if electrode_name == "MicroElectrode":
        #     desired.add("fillet")

        assert desired == set(faces)

    def check_electrode_volume(self, electrode, electrode_name):
        """Check volume of the entire electrode."""
        desired = self._calculate_electrode_volume(electrode, electrode_name)
        actual = electrode.geometry.mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def check_contacts_volume(self, electrode, electrode_name):
        """Check volume of all the contacts."""
        desired = self._calculate_contacts_volume(electrode, electrode_name)
        actual = electrode._contacts().mass
        tolerance = 1e-5

        np.testing.assert_allclose(actual, desired, atol=tolerance)

    def _calculate_electrode_volume(self, electrode, electrode_name):
        if electrode_name == "MicroProbesRodentElectrode":
            total_length = electrode._parameters.total_length
            lead_radius = electrode._parameters.lead_radius
            contact_radius = electrode._parameters.contact_radius
            height = total_length - contact_radius

            return (height * lead_radius**2 * np.pi) + (
                4 / 3 * np.pi * contact_radius**3 * 0.5
            )

        elif electrode_name == "MicroProbesSNEX100":
            radius_1 = electrode._parameters.core_tubing_diameter * 0.5
            height_1 = electrode._parameters.core_tubing_length

            distance_2 = (
                electrode._parameters.core_electrode_length
                + electrode._parameters.core_tubing_length
                + electrode._parameters.outer_electrode_length
            )

            radius_2 = electrode._parameters.outer_tubing_diameter * 0.5
            height_2 = electrode._parameters.total_length - distance_2

            body1 = np.pi * radius_1**2 * height_1
            body2 = np.pi * radius_2**2 * height_2

            outer_electrode_radius = (
                electrode._parameters.outer_electrode_diameter * 0.5
            )
            outer_electrode_length = electrode._parameters.outer_electrode_length

            core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
            core_electrode_length = electrode._parameters.core_electrode_length

            contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
            contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
                np.pi
                * core_electrode_radius**2
                * (core_electrode_length - core_electrode_radius)
            )

            return body1 + body2 + contact_1 + contact_2

        elif "NeuroNexus" in electrode_name:
            total_length = electrode._parameters.total_length
            tip_length = electrode._parameters.tip_length
            contact_spacing = electrode._parameters.contact_spacing
            first_contact_position = electrode._parameters.get_center_first_contact()
            thickness = electrode._parameters.shank_thickness
            max_width = electrode._parameters.max_width
            min_width = electrode._parameters.min_width
            n_contacts = electrode._n_contacts
            # length of straight part until first contact, where width
            # is max_width
            straight_part_length = (
                total_length - tip_length - (n_contacts - 1) * contact_spacing
            )
            shaft_volume = straight_part_length * max_width * thickness

            # surface area until tip, shaped like a trapezoid
            middle_area = (
                0.5 * (min_width + max_width) * (n_contacts - 1) * contact_spacing
            )
            middle_volume = middle_area * thickness

            # last tip part
            tip_area = 0.5 * min_width * first_contact_position
            tip_volume = tip_area * thickness
            return shaft_volume + middle_volume + tip_volume

        else:
            contact_length = electrode._parameters.contact_length
            lead_radius = electrode._parameters.lead_diameter * 0.5
            tip_length = electrode._parameters.tip_length
            total_length = electrode._parameters.total_length
            height = total_length - tip_length

            if electrode_name == "MicroElectrode":
                tip_radius = electrode._parameters.tip_diameter * 0.5
                filet_val = 0.001142187023875807

                return (
                    (contact_length * tip_radius**2 * np.pi)
                    + (height * lead_radius**2 * np.pi)
                    - filet_val
                )

            else:
                return (np.pi * lead_radius**2 * height) + (
                    4 / 3 * np.pi * lead_radius**3 * 0.5
                )

    def _calculate_contacts_volume(self, electrode, electrode_name):
        if electrode_name == "MicroProbesRodentElectrode":
            contact_radius = electrode._parameters.contact_radius
            return 4 / 3 * np.pi * contact_radius**3 * 0.5
        elif electrode_name == "MicroProbesSNEX100":
            outer_electrode_radius = (
                electrode._parameters.outer_electrode_diameter * 0.5
            )
            outer_electrode_length = electrode._parameters.outer_electrode_length

            core_electrode_radius = electrode._parameters.core_electrode_diameter * 0.5
            core_electrode_length = electrode._parameters.core_electrode_length

            contact_1 = np.pi * outer_electrode_radius**2 * outer_electrode_length
            contact_2 = (4 / 3 * np.pi * core_electrode_radius**3 * 0.5) + (
                np.pi
                * core_electrode_radius**2
                * (core_electrode_length - core_electrode_radius)
            )

            return contact_1 + contact_2
        else:
            contact_length = electrode._parameters.contact_length
            lead_radius = electrode._parameters.lead_diameter * 0.5
            tip_length = electrode._parameters.tip_length
            n_contacts = electrode._n_contacts

            if electrode_name == "BostonScientificCartesiaHX":
                return (np.pi * lead_radius**2 * contact_length) * 4 + (
                    np.pi * lead_radius**2 * contact_length * 90 / 360
                ) * 12
            elif electrode_name == "BostonScientificCartesiaX":
                return (np.pi * lead_radius**2 * contact_length) + (
                    np.pi * lead_radius**2 * contact_length * 90 / 360
                ) * 15
            elif electrode_name == "BostonScientificVerciseDirected":
                C1_height = tip_length - lead_radius
                C1_volume = (4 / 3 * lead_radius**3 * np.pi * 0.5) + (
                    C1_height * lead_radius**2 * np.pi
                )
                return (
                    C1_volume
                    + (np.pi * lead_radius**2 * contact_length)
                    + (np.pi * lead_radius**2 * contact_length * 90 / 360) * 6
                )

            elif (
                electrode_name == "BostonScientificVercise"
                or electrode_name in self.Medtronic
                or electrode_name in self.NeuroPace
                or electrode_name in self.PINSMedical
                or electrode_name in self.SceneRay1242
            ):
                C1_height = tip_length - lead_radius
                C1_volume = (4 / 3 * lead_radius**3 * np.pi * 0.5) + (
                    C1_height * lead_radius**2 * np.pi
                )

                return (contact_length * lead_radius**2 * np.pi) * n_contacts
            elif (
                electrode_name in self.AbbottStJudeDirected
                or electrode_name in self.MedtronicSenSight
            ):
                return (np.pi * lead_radius**2 * contact_length * 2) + (
                    np.pi * lead_radius**2 * contact_length * 90 / 360
                ) * 6
            elif electrode_name == "MicroElectrode":
                tip_radius = electrode._parameters.tip_diameter * 0.5
                filet_val = 0.00114218702387580
                return (contact_length * tip_radius**2 * np.pi) - filet_val
            else:
                if electrode_name in self.Dixi or electrode_name in self.PMTsEEG2102:
                    C1_height = contact_length - lead_radius
                else:
                    C1_height = tip_length - lead_radius
                C1_volume = (4 / 3 * lead_radius**3 * np.pi * 0.5) + (
                    C1_height * lead_radius**2 * np.pi
                )

                return (contact_length * lead_radius**2 * np.pi) * (
                    n_contacts - 1
                ) + C1_volume
