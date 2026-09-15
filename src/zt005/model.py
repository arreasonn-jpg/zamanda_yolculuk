from dataclasses import dataclass

@dataclass(frozen=True)
class ActiveCylinder:
    radius_m: float = 1.0
    height_m: float = 2.5

@dataclass(frozen=True)
class Backpack:
    height_m: float = 0.45
    width_m: float = 0.32
    depth_m: float = 0.15
    center_y_m: float = -1.075

@dataclass(frozen=True)
class Drive:
    frequency_hz: float = 1.0e6
    peak_current_a: float = 100.0
