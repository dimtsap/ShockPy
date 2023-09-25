from abc import ABC, abstractmethod

import numpy as np
from scipy.interpolate import splprep, splev
from shapely.geometry import LineString

from ImpedancePy.shock_wave_compression.material_states.Hugoniot import Hugoniot
from ImpedancePy.shock_wave_compression.material_states.Intersection import Intersection
from ImpedancePy.shock_wave_compression.material_states.Isentrope import Isentrope
from ImpedancePy.shock_wave_compression.material_states.isentrope_calculators.IntegratedIsentrope import IntegratedIsentrope


class Material(ABC):
    isentrope_calculator = IntegratedIsentrope()

    def __init__(self, gamma_eff: float,
                 initial_density: float,
                 released=True,
                 is_stochastic=False):
        """

        :param gamma_eff:
        :param initial_density:
        :param released:
        :param is_stochastic:
        """
        self.is_stochastic = is_stochastic
        self.initial_density = initial_density
        self.initial_volume = 1 / initial_density
        self.Gamma_eff = gamma_eff
        self.isentropes = []
        self.hugoniots_list = []
        self.nominal_hugoniot: Hugoniot = None
        self.released = released

    def initialize_isentropes_at_pressure(self, shock_pressure: float):
        point_hugoniot_tuples = self.release_isentropes_at_pressure(shock_pressure)
        release_isentropes = [Material.isentrope_calculator.calculate_isentrope(intersection, hugoniot, self.Gamma_eff,
                                                                                self.initial_volume, self.released)
                              for (intersection, hugoniot) in point_hugoniot_tuples]
        return release_isentropes

    def calculate_intersection(self, hugoniot, isentrope):
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
            intersection = hugoniot.find_hugoniot_point_at_pressure(pressure, self.initial_density)
            release_isentrope = self.isentrope_calculator.calculate_isentrope(hugoniot, intersection, self.Gamma_eff,
                                                                              self.initial_volume, self.released)
            isentropes.append(release_isentrope)
        return isentropes

    def _find_hugoniot_point_at_shock_velocity(self, hugoniot: Hugoniot, shock_velocity: float):
        self.Gamma_eff = 0.619 * (1 - np.exp(-0.0882 * (shock_velocity - 12.0922) ** 1.5))
        particle_velocity = hugoniot.interpolate_shock_velocity(shock_velocity)
        current_density = self.initial_density * shock_velocity / (
                shock_velocity - particle_velocity)
        current_volume = 1 / current_density
        compression_ratio = self.initial_volume / current_volume
        pressure = hugoniot.interpolate_particle_velocity(particle_velocity)
        return Intersection(pressure=pressure, particle_velocity=particle_velocity,
                            shock_velocity=shock_velocity, volume=current_volume,
                            compression_ratio=compression_ratio,
                            hugoniot=hugoniot)

    def hugoniots_intersection_with_isentrope(self, isentrope: Isentrope):
        new_isentropes = []
        for hugoniot in self.hugoniots_list:
            intersection = self.calculate_intersection(hugoniot, isentrope)
            release_isentrope = self.isentrope_calculator.calculate_isentrope(hugoniot, intersection, self.Gamma_eff,
                                                                              self.initial_volume, self.released)
            new_isentropes.append(release_isentrope)
        return new_isentropes

    @abstractmethod
    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        pass

    def calculate_hugoniot(self, parameters):
        hugoniot_particle_velocity = np.linspace(0, 30, num=1000)
        hugoniot_shock_velocity = self.analytical_shock_velocity_equation(parameters, hugoniot_particle_velocity)
        hugoniot_pressure = self.initial_density * hugoniot_shock_velocity * hugoniot_particle_velocity
        hugoniot_volume = self.initial_volume * (
                hugoniot_shock_velocity - hugoniot_particle_velocity) / hugoniot_shock_velocity
        return Hugoniot(pressures=hugoniot_pressure,
                        particle_velocities=hugoniot_particle_velocity,
                        shock_velocities=hugoniot_shock_velocity,
                        volumes=hugoniot_volume)

    def intersections_at_pressure(self, pressure):
        intersections = [hugoniot.find_hugoniot_point_at_pressure(pressure, self.initial_density)
                         for hugoniot in self.hugoniots_list]
        return intersections

    def find_previous_material_intersections_from_current(self, current_intersection: Intersection, next_material):
        intersections = []
        for hugoniot in self.hugoniots_list:
            intersection = self.isentrope_calculator \
                .find_previous_material_isentrope(previous_material_hugoniot=hugoniot,
                                                  current_material_intersection=current_intersection,
                                                  previous_material=self,
                                                  next_material=next_material)
            intersections.append(intersection)
        return intersections
