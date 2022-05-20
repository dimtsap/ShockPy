from dataclasses import dataclass
from shock_wave_compression.material_states.Hugoniot import Hugoniot


@dataclass
class Intersection:
    pressure: float
    particle_velocity: float
    shock_velocity: float
    volume: float
    compression_ratio: float
    hugoniot: Hugoniot = None
