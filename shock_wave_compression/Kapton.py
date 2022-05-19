import numpy as np
from scipy import interpolate

from shock_wave_compression.Material import Material


class Kapton(Material):
    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 1.33):  # Source https://dielectricmfg.com/knowledge-base/kapton/
        super().__init__(1.1, initial_density)
        self.hugoniot_P_up = np.loadtxt('kapton_hugoniot.txt')
        self.hugoniot_pressure = self.hugoniot_P_up[:, 0]
        self.hugoniot_particle_velocity = self.hugoniot_P_up[:, 1]
        self.hugoniot_shock_velocity = self.hugoniot_pressure / (self.initial_density * self.hugoniot_particle_velocity)
        self.hugoniot_volume = self.initial_volume * (
                self.hugoniot_shock_velocity - self.hugoniot_particle_velocity) / self.hugoniot_shock_velocity
        self._hugoniot_interpolator = interpolate.interp1d(self.hugoniot_P_up[:, 0],
                                                           self.hugoniot_P_up[:, 1])

    def calculate_nominal_hugoniot(self):
        return (self.hugoniot_pressure, self.hugoniot_particle_velocity,
                self.hugoniot_shock_velocity, self.hugoniot_volume)

    def calculate_isentrope(self, input_pressure,
                            input_volume,
                            input_particle_velocity,
                            compression_ratio):
        release_pressures, release_up = super().calculate_isentrope(input_pressure, input_volume,
                                                                    input_particle_velocity, compression_ratio)
        # tries to extrapolate Kapton isentrope
        f = interpolate.interp1d(release_up[0][1:], release_pressures[0][1:], kind='slinear',
                                 fill_value='extrapolate')
        release_kapton_up = np.linspace(2, 15)
        release_kapton_p = f(release_kapton_up)

        return [release_kapton_p], [release_kapton_up]

    def find_hugoniot_point_at_pressure(self, pressure: float):
        particle_velocity = self._hugoniot_interpolator(pressure)
        shock_velocity = pressure / (self.initial_density * particle_velocity)  # P=rho0*Us*up
        current_density = self.initial_density * shock_velocity / (
                shock_velocity - particle_velocity)  # rho0*Us=rho1*(Us-up)
        current_volume = 1 / current_density
        compression_ratio = self.initial_volume / current_volume
        return particle_velocity, current_volume, compression_ratio
