# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Package to prepare geometries for VCM simulation."""

from .bounding_box import BoundingBox
from .contacts import Contact, Contacts
from .model_geometry import BrainGeometry, ModelGeometry

__all__ = ("BoundingBox", "BrainGeometry", "Contact", "Contacts", "ModelGeometry")
