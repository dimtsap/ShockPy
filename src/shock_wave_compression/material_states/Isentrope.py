from dataclasses import dataclass
import numpy as np


@dataclass
class Isentrope:
    pressures: np.ndarray
    particle_velocities: np.ndarray
    intersection: object
