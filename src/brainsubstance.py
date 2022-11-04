from enum import IntEnum, unique


@unique
class Material(IntEnum):
    UNKNOWN = 0
    CEREBRO_SPINAL_FLUID = 1
    WHITE_MATTER = 2
    GRAY_MATTER = 3
