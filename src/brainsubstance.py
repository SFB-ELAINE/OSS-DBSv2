from enum import IntEnum, unique


@unique
class BrainSubstance(IntEnum):
    UNKNOWN = 0
    CEREBROSPINAL_FLUID = 1
    WHITE_MATTER = 2
    GRAY_MATTER = 3
