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

    def interpolate_pressure(self, pressure: float) -> float:
        hugoniot_interpolator = interpolate.interp1d(self.pressures, self.particle_velocities)
        return hugoniot_interpolator(pressure)

    def interpolate_shock_velocity(self, shock_velocity: float) -> float:
        hugoniot_interpolator = interpolate.interp1d(self.shock_velocities, self.particle_velocities)
        return hugoniot_interpolator(shock_velocity)

    def interpolate_particle_velocity(self, particle_velocity: float) -> float:
        hugoniot_interpolator = interpolate.interp1d(self.particle_velocities, self.pressures)
        return hugoniot_interpolator(particle_velocity)