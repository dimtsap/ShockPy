from dataclasses import dataclass
import numpy as np

from shock_wave_compression.material_states.Intersection import Intersection


@dataclass
class Isentrope:
    pressures: np.ndarray
    particle_velocities: np.ndarray
    intersection: Intersection
