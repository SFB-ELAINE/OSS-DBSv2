from enum import IntEnum, unique


@unique
class Material(IntEnum):
    """Enumeration for brain substances."""
    UNKNOWN = 0
    GRAY_MATTER = 1
    WHITE_MATTER = 2
    CSF = 3
    BLOOD = 4
