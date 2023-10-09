from abc import ABC, abstractmethod

import numpy as np
from scipy.interpolate import splprep, splev
from shapely.geometry import LineString

from ShockPy.shock_wave_compression.material_states.Hugoniot import Hugoniot
from ShockPy.shock_wave_compression.material_states.Intersection import Intersection
from ShockPy.shock_wave_compression.material_states.Isentrope import Isentrope
from ShockPy.shock_wave_compression.material_states.isentrope_calculators.IntegratedIsentrope import IntegratedIsentrope


class Material(ABC):
    isentrope_calculator = IntegratedIsentrope()
    """Static attribute common for all materials throughout an experiment. Defines the way in which the code 
    approximates the release isentrope of a material, with two available alternatives, either 
    :class:`.ReflectedHugoniot` or :class:`.IntegratedIsentrope`"""

    def __init__(self, gamma_eff: float,
                 initial_density: float,
                 released=True,
                 is_stochastic=False):
        """
        Initializer of the Material baseclass.

        :param gamma_eff: :math:`\Gamma_{eff}` is the updated Mie-Grüneisen (MG) parameter :math:`\Gamma`
         for the linear reference of the Mie-Grüneisen (MGLR) model, that assumes the release path can be accurately
         reproduced with constant :math:`\Gamma` along nearly its entirety. More details on the process of calculating
         the :math:`\Gamma_{eff}` given experimental data can be found in :cite:`Knudson`.
        :param initial_density: Ambient density of the material.
        :param released: Boolean parameter indicating whether the reflected shock produced at the interface of the
         current material with the next is a rarefaction or a reshock.
          * If :any:`True`, then the algorithms assumes that the material undergoes release and compute the release part of the isentrope.
          * If :any:`False`, then the algorithm assumes the material will be reshocked and computed the respective part of the isentrope.
         Default value: :any:`True`.
        :param is_stochastic: Boolean that defines if the Hugoniot used for this material is the deterministic one or
         multiple Hugoniots that are produced as a result of the Bayesian Inference framework will be used.
          * If :any:`True`, then samples of the Bayesian Inference will be used to generate the :code:`hugoniots_list` attributed.
          * If :any:`False`, only the nominal Hugoniot populates the list, and is derived from the least-square fitting of its analytical equation to the experimental data.
         Default value: :any:`False`.
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
        """
        Auxiliary function:
        Given an input shock pressure, it finds the shocked Hugoniot state of the material and subsequently the
        release Isentrope that goes through the specific shocked state. In case of an uncertain Hugoniot, a single
        isentrope for each Hugoniot is calculated.

        :param shock_pressure: Input laser drive (GPa)
        :return: A list containing the Isentropes of the materia at a specific input pressure.
        """
        point_hugoniot_tuples = self.release_isentropes_at_pressure(shock_pressure)
        release_isentropes = [Material.isentrope_calculator.calculate_isentrope(intersection, hugoniot, self.Gamma_eff,
                                                                                self.initial_volume, self.released)
                              for (intersection, hugoniot) in point_hugoniot_tuples]
        return release_isentropes

    def calculate_intersection(self, hugoniot: Hugoniot, isentrope: Isentrope):
        """
        Given a :class:`.Hugoniot` and an :class:`.Isentrope` this function calculates their intersection in the
        :math:`(u_p - P)` space and then calculates all other quantities (e.g. :math:`U_s`) using the Rankine Hugoniot
        equations.

        :param hugoniot: A Hugoniot object containing all required info about the next material.
        :param isentrope: An Isentrope object of the current material
        :return: An :class:`.Intersection` containing the common state of the Hugoniot and Intersection
        """
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
        """
        Auxiliary function for the backward propagation of the experiment. Given the current material shocked state in
        the form of an intersection, as well as a reference to the next material, tries to find the Isentrope of the
        current material that goes through the next material shocked state.

        :param current_intersection: An :class:`.Intersection` object containing the next material shocked state.
        :param next_material: A reference to the next material which gives access to its Hugoniots and material
         properties
        :return: A list containing one (if the current material is deterministic) or multiple shocked states
         (if the material is uncertain) in the form of :class:`.Intersection` objects.
        """
        intersections = []
        for hugoniot in self.hugoniots_list:
            intersection = self.isentrope_calculator \
                .find_previous_material_isentrope(previous_material_hugoniot=hugoniot,
                                                  current_material_intersection=current_intersection,
                                                  previous_material=self,
                                                  next_material=next_material)
            intersections.append(intersection)
        return intersections
