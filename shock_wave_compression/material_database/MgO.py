import numpy as np
from scipy import interpolate
from shock_wave_compression.material_database.baseclass.Material import Material
from shock_wave_compression.material_states.Hugoniot import Hugoniot


class MgO(Material):
    def release_isentropes_at_pressure(self, pressure: float):
        pass

    def calculate_stochastic_hugoniot(self) -> list[Hugoniot]:
        pass

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 3.58):
        super().__init__(gamma_eff, initial_density)
        self.hugoniot_particle_velocity = np.linspace(0, 20, num=1000)
        self.hugoniot_shock_velocity = 6.6161 + 1.4111 * self.hugoniot_particle_velocity - 0.016277 * np.square(
            self.hugoniot_particle_velocity)
        self.hugoniot_pressure = self.initial_density * self.hugoniot_shock_velocity * self.hugoniot_particle_velocity
        self.hugoniot_volume = self.initial_volume * (
                self.hugoniot_shock_velocity - self.hugoniot_particle_velocity) / self.hugoniot_shock_velocity
        self._hugoniot_interpolator = interpolate.interp1d(self.hugoniot_pressure,
                                                           self.hugoniot_particle_velocity)

    def calculate_nominal_hugoniot(self):
        return (self.hugoniot_pressure, self.hugoniot_particle_velocity,
                self.hugoniot_shock_velocity, self.hugoniot_volume)

    def find_hugoniot_point_at_pressure(self, pressure: float):
        particle_velocity = self._hugoniot_interpolator(pressure)
        shock_velocity = pressure / (self.initial_density * particle_velocity)  # P=rho0*Us*up
        current_density = self.initial_density * shock_velocity / (
                shock_velocity - particle_velocity)  # rho0*Us=rho1*(Us-up)
        current_volume = 1 / current_density
        compression_ratio = self.initial_volume / current_volume
        return particle_velocity, current_volume, compression_ratio
