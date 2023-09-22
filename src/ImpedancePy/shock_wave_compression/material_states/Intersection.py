from dataclasses import dataclass

from ImpedancePy.shock_wave_compression.material_states.Hugoniot import Hugoniot
from ImpedancePy.shock_wave_compression.material_states.Isentrope import Isentrope


@dataclass
class Intersection:
    pressure: float
    particle_velocity: float
    shock_velocity: float
    volume: float
    compression_ratio: float
    hugoniot: Hugoniot = None
    isentrope: Isentrope = None
