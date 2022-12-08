from abc import ABC, abstractmethod

import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import splprep, splev
from shapely.geometry import LineString

from shock_wave_compression.material_states.Hugoniot import Hugoniot
from shock_wave_compression.material_states.Intersection import Intersection
from shock_wave_compression.material_states.Isentrope import Isentrope


class Material(ABC):

    def __init__(self, gamma_eff: float,
                 initial_density: float,
                 is_stochastic=False):
        self.is_stochastic = is_stochastic
        self.initial_density = initial_density
        self.initial_volume = 1 / initial_density
        self.Gamma_eff = gamma_eff
        self.isentropes = []
        self.hugoniots_list = []
        self.nominal_hugoniot: Hugoniot = None

    def initialize_isentropes_at_pressure(self, shock_pressure: float):
        point_hugoniot_tuples = self.release_isentropes_at_pressure(shock_pressure)
        release_isentropes = [self.calculate_isentrope(intersection, hugoniot)
                              for (intersection, hugoniot) in point_hugoniot_tuples]
        return release_isentropes

    def calculate_isentrope(self, hugoniot, intersection):
        pressuresList = hugoniot.pressures.tolist()
        volumesList = hugoniot.volumes.tolist()
        release_pressures = []
        release_particle_velocities = []
        for (volume_Hugoniot, pressure_Hugoniot) in zip([volumesList], [pressuresList]):
            # V = np.linspace(0.01, intersection.volume, num=1000)
            V = np.linspace(intersection.volume, max(volume_Hugoniot), num=1000)

            S = 1.31867  # What is this parameter.? Is it material dependent.?
            C01 = intersection.pressure / (self.initial_density * intersection.particle_velocity) \
                  - S * intersection.particle_velocity

            reference_hugoniot_pressure = self.initial_density * C01 ** 2 * (self.initial_volume / V - 1) \
                                          * (self.initial_volume / V) / np.square(S - (S - 1) * self.initial_volume / V)

            energy_integral = []
            for i in range(len(V)):
                delta_energy_integral_function = lambda x: (x / intersection.volume) ** self.Gamma_eff * \
                                                           reference_hugoniot_pressure[i] \
                                                           * (1 - self.Gamma_eff / 2 * (self.initial_volume / x - 1))
                aux = quad(delta_energy_integral_function, intersection.volume, V[i])
                energy_integral.append(aux[0])

            Es_E0 = [intersection.pressure * self.initial_volume / 2 * (
                    intersection.compression_ratio - 1) / intersection.compression_ratio * (
                             intersection.volume / x) ** self.Gamma_eff - (
                             intersection.volume / x) ** self.Gamma_eff * y
                     for (x, y) in zip(V, energy_integral)]

            release_pressure = [PH * (1 - self.Gamma_eff / 2 * (self.initial_volume / v - 1)) + self.Gamma_eff / v * dE
                                for (PH, v, dE) in zip(reference_hugoniot_pressure, V, Es_E0)]

            release_particle_velocity = [float(intersection.particle_velocity)]

            ups = list(np.real(
                intersection.particle_velocity + cumulative_trapezoid(np.sqrt(-np.diff(release_pressure) / np.diff(V)),
                                                                      V[1:],
                                                                      initial=0)))
            release_particle_velocity.extend(ups)

            release_pressures.append(np.array(release_pressure))
            release_particle_velocities.append(np.array(release_particle_velocity))

        return Isentrope(pressures=np.atleast_1d(release_pressures),
                         particle_velocities=np.atleast_1d(release_particle_velocities),
                         intersection=intersection)

    def calculate_intersection(self, hugoniot, isentrope):
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(hugoniot.particle_velocities, hugoniot.pressures)
        # plt.plot(np.squeeze(isentrope.particle_velocities), np.squeeze(isentrope.pressures))
        # plt.show()
        line_1 = LineString(np.column_stack((hugoniot.particle_velocities, hugoniot.pressures)))
        line_2 = LineString(np.column_stack((np.squeeze(isentrope.particle_velocities),
                                             np.squeeze(isentrope.pressures))))
        intersection = line_1.intersection(line_2)
        pressure = intersection.xy[1][0]
        particle_velocity = intersection.xy[0][0]
        shock_velocity = pressure / (self.initial_density * particle_velocity)
        current_density = self.initial_density * shock_velocity / (
                shock_velocity - particle_velocity)  # rho0*Us=rho1*(Us-up)
        current_volume = 1 / current_density
        compression_ratio = self.initial_volume / current_volume
        return Intersection(pressure=pressure,
                            particle_velocity=particle_velocity,
                            shock_velocity=shock_velocity,
                            volume=current_volume,
                            compression_ratio=compression_ratio,
                            hugoniot=hugoniot,
                            isentrope=isentrope)

    @staticmethod
    def _upsample_coords(coord_list):
        # s is smoothness, set to zero
        # k is degree of the spline. setting to 1 for linear spline
        tck, u = splprep(coord_list, k=1, s=0.0)
        upsampled_coords = splev(np.linspace(0, 1, 100), tck)
        return upsampled_coords

    def release_isentropes_at_pressure(self, pressure: float):
        isentropes = []
        for hugoniot in self.hugoniots_list:
            intersection = self._find_hugoniot_point_at_pressure(hugoniot, pressure)
            release_isentrope = self.calculate_isentrope(hugoniot, intersection)
            isentropes.append(release_isentrope)
        return isentropes

    def _find_hugoniot_point_at_pressure(self, hugoniot: Hugoniot, pressure: float):
        particle_velocity = hugoniot.interpolate(pressure)
        shock_velocity = pressure / (self.initial_density * particle_velocity)  # P=rho0*Us*up
        current_density = self.initial_density * shock_velocity / (
                shock_velocity - particle_velocity)  # rho0*Us=rho1*(Us-up)
        current_volume = 1 / current_density
        compression_ratio = self.initial_volume / current_volume
        return Intersection(pressure=pressure, particle_velocity=particle_velocity,
                            shock_velocity=shock_velocity, volume=current_volume,
                            compression_ratio=compression_ratio,
                            hugoniot=hugoniot)

    def hugoniots_intersection_with_isentrope(self, isentrope: Isentrope):
        new_isentropes = []
        for hugoniot in self.hugoniots_list:
            intersection = self.calculate_intersection(hugoniot, isentrope)
            release_isentrope = self.calculate_isentrope(hugoniot, intersection)
            new_isentropes.append(release_isentrope)
        return new_isentropes

    @abstractmethod
    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        pass

    def calculate_hugoniot(self, parameters):
        hugoniot_particle_velocity = np.linspace(0, 20, num=1000)
        hugoniot_shock_velocity = self.analytical_shock_velocity_equation(parameters, hugoniot_particle_velocity)
        hugoniot_pressure = self.initial_density * hugoniot_shock_velocity * hugoniot_particle_velocity
        hugoniot_volume = self.initial_volume * (
                hugoniot_shock_velocity - hugoniot_particle_velocity) / hugoniot_shock_velocity
        return Hugoniot(pressures=hugoniot_pressure,
                        particle_velocities=hugoniot_particle_velocity,
                        shock_velocities=hugoniot_shock_velocity,
                        volumes=hugoniot_volume)
