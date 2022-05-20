from dataclasses import dataclass
import numpy as np
from scipy.interpolate import interpolate


@dataclass
class Hugoniot:
    pressure: np.ndarray
    particle_velocity: np.ndarray
    shock_velocity: np.ndarray
    volume: np.ndarray
    id: int = None

    def interpolate(self, pressure: float) -> float:
        hugoniot_interpolator = interpolate.interp1d(self.pressure, self.particle_velocity)
        return hugoniot_interpolator(pressure)
