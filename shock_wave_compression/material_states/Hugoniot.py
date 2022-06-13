from dataclasses import dataclass
import numpy as np
from scipy.interpolate import interpolate


@dataclass
class Hugoniot:
    pressures: np.ndarray
    particle_velocities: np.ndarray
    shock_velocities: np.ndarray
    volumes: np.ndarray
    id: int = None

    def interpolate(self, pressure: float) -> float:
        hugoniot_interpolator = interpolate.interp1d(self.pressures, self.particle_velocities)
        return hugoniot_interpolator(pressure)
