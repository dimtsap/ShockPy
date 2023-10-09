from dataclasses import dataclass
import numpy as np
from scipy.interpolate import interpolate


@dataclass
class Hugoniot:
    """
    Dataclass containing all info (e.g. :math:`P, u_p, U_s, V`) of a Hugoniot curve.

    :param pressures: A numpy array containing the pressures :math:`(P)` for all shocked states of a Hugoniot.
    :param particle_velocities: A numpy array containing the particle velocities :math:`(u_p)` for all shocked states
     of a Hugoniot.
    :param shock_velocities: A numpy array containing the shock velocities :math:`(U_s)` for all shocked states
     of a Hugoniot.
    :param volumes: A numpy array containing the pressures for all shocked states of a Hugoniot.
    :param id: An integer identifier.
    """
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

    def reflected_hugoniot(self, initial_density, intersection_particle_velocity=0, delta_particle_velocity=0):
        """
        Calculates and returns the reflection of a Hugoniot given either a vertical axis or a shift value.

        :param initial_density: Ambient density of the material
        :param intersection_particle_velocity: Particle velocity value that defines the vertical reflection axis of the
         Hugoniot.
        :param delta_particle_velocity: Value of particle velocity to make horizontal shifts to the reflected Hugoniot.
        """
        reflected_hugoniot_particle_velocities = intersection_particle_velocity + \
                                                 (intersection_particle_velocity - self.particle_velocities) \
                                                 + delta_particle_velocity
        initial_volume = 1 / initial_density
        us = self.pressures / (initial_density * reflected_hugoniot_particle_velocities)
        v = initial_volume * (us - reflected_hugoniot_particle_velocities) / us

        return Hugoniot(
            pressures=self.pressures.__copy__(),
            particle_velocities=reflected_hugoniot_particle_velocities,
            shock_velocities=us,
            volumes=v
        )

    def find_hugoniot_point_at_pressure(self, pressure: float, initial_density: float):
        initial_volume = 1 / initial_density
        particle_velocity = self.interpolate_pressure(pressure)
        shock_velocity = pressure / (initial_density * particle_velocity)  # P=rho0*Us*up
        current_density = initial_density * shock_velocity / (
                shock_velocity - particle_velocity)  # rho0*Us=rho1*(Us-up)
        current_volume = 1 / current_density
        compression_ratio = initial_volume / current_volume
        from src.ShockPy.shock_wave_compression.material_states.Intersection import Intersection
        return Intersection(pressure=pressure, particle_velocity=particle_velocity,
                            shock_velocity=shock_velocity, volume=current_volume,
                            compression_ratio=compression_ratio,
                            hugoniot=self)
