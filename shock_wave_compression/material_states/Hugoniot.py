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

    def reflected_hugoniot(self, initial_density, intersection_particle_velocity=0, delta_particle_velocity=0):
        reflected_hugoniot_particle_velocities = intersection_particle_velocity + \
                                                 (intersection_particle_velocity - self.particle_velocities)\
                                                 +delta_particle_velocity
        initial_volume = 1 / initial_density
        us = self.pressures / (initial_density * reflected_hugoniot_particle_velocities)
        v = initial_volume * (us - reflected_hugoniot_particle_velocities) / us

        return Hugoniot(
            pressures=self.pressures.__copy__(),
            particle_velocities=reflected_hugoniot_particle_velocities,
            shock_velocities=us,
            volumes=v
        )
