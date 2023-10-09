import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import interpolate

from ShockPy.shock_wave_compression.material_states.Hugoniot import Hugoniot
from ShockPy.shock_wave_compression.material_states.Intersection import Intersection
from ShockPy.shock_wave_compression.material_states.Isentrope import Isentrope


class IntegratedIsentrope:

    def calculate_isentrope(self, hugoniot, intersection, Gamma_eff, initial_volume, released=True):
        """
        An numerical integration approach to calculating the release isentrope of a material based on :cite:`Knudson`
        method for Quartz.

        :param hugoniot: Hugoniot of the current material
        :param intersection: An :class:`.Intersection` object containing the information about the current material
         shocked state.
        :param Gamma_eff: The :math:`\Gamma_{eff}` Gruneisen parameter of the current material
        :param initial_volume: The ambient volume of the current material
        :param released: released: Boolean parameter indicating whether the reflected shock produced at the interface
         of the current material with the next is a rarefaction or a reshock.
        """
        pressuresList = hugoniot.pressures.tolist()
        volumesList = hugoniot.volumes.tolist()
        release_pressures = []
        release_particle_velocities = []
        V = np.linspace(intersection.volume, max(volumesList), num=1000) if released \
            else np.linspace(intersection.volume, 0.001, num=1000)
        hugoniot_interpolator = interpolate.interp1d(np.flip(np.array(volumesList)), np.flip(np.array(pressuresList)),
                                                     kind='cubic',
                                                     fill_value="extrapolate")
        reference_hugoniot_pressure = hugoniot_interpolator(V)

        energy_integral = []
        for i in range(len(V)):
            delta_energy_integral_function = lambda x: (x / intersection.volume) ** Gamma_eff * \
                                                       reference_hugoniot_pressure[i] \
                                                       * (1 - Gamma_eff / 2 * (initial_volume / x - 1))
            aux = quad(delta_energy_integral_function, intersection.volume, V[i])
            energy_integral.append(aux[0])

        Es_E0 = [intersection.pressure * initial_volume / 2 * (
                intersection.compression_ratio - 1) / intersection.compression_ratio * (
                         intersection.volume / x) ** Gamma_eff - (
                         intersection.volume / x) ** Gamma_eff * y
                 for (x, y) in zip(V, energy_integral)]

        release_pressure = [PH * (1 - Gamma_eff / 2 * (initial_volume / v - 1)) + Gamma_eff / v * dE
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

    def find_previous_material_isentrope(self,
                                         previous_material_hugoniot: Hugoniot,
                                         current_material_intersection: Intersection,
                                         previous_material, next_material):
        """
        This methods evaluates the Isentrope of the previous material, that goes through the current material shocked
        state. The two material Hugoniots and the current material shocked state are required. This functions aids in
        the backward propagation of the experiment from the window to the input laser drive. SIne the evaluation of
        the isentrope in this case is by numerical integration, a secant method is utilized to find the appropriate
        isentrope of the previous material.


        :param previous_material_hugoniot: The Hugoniot of the previous material. In both Isentrope Calculators it is
         useful for obtaining the required isentrope and thus propagating the experiment backwards.
        :param current_material_intersection: An :class:`Intersection` object that contains the information of the
         current material shocked state.
        :param previous_material: An object reference to the previous material and its properties
        :param next_material: An object reference to the current material and its properties:
        """
        previous_material_density = previous_material.initial_density
        min_intersection = previous_material_hugoniot \
            .find_hugoniot_point_at_pressure(100, previous_material_density)
        max_intersection = previous_material_hugoniot \
            .find_hugoniot_point_at_pressure(700, previous_material_density)

        intersection, iterations = self.__secant_method(min_intersection, max_intersection,
                                                        previous_material_hugoniot, current_material_intersection,
                                                        previous_material, next_material)
        return intersection

    def __secant_method(self, intersection0, intersection1,
                        previous_material_hugoniot,
                        next_material_intersection,
                        previous_material, next_material,
                        tol=1e-5, n=0):
        # increment counter
        n += 1

        next_material_hugoniot = next_material_intersection.hugoniot
        # calculate function values at endpoints
        isentrope0 = self.calculate_isentrope(hugoniot=previous_material_hugoniot,
                                              intersection=intersection0,
                                              Gamma_eff=previous_material.Gamma_eff,
                                              initial_volume=previous_material.initial_volume,
                                              released=previous_material.released)
        goal_intersection0 = next_material.calculate_intersection(hugoniot=next_material_hugoniot,
                                                                  isentrope=isentrope0)

        isentrope1 = self.calculate_isentrope(hugoniot=previous_material_hugoniot,
                                              intersection=intersection1,
                                              Gamma_eff=previous_material.Gamma_eff,
                                              initial_volume=previous_material.initial_volume,
                                              released=previous_material.released)

        goal_intersection1 = next_material.calculate_intersection(hugoniot=next_material_hugoniot,
                                                                  isentrope=isentrope1)

        # calculate next root approximation
        x1 = intersection1.pressure
        y1 = goal_intersection1.pressure - next_material_intersection.pressure
        x0 = intersection0.pressure
        y0 = goal_intersection0.pressure - next_material_intersection.pressure

        xn = x1 - y1 * ((x1 - x0) / (y1 - y0))

        new_intersection = previous_material_hugoniot.find_hugoniot_point_at_pressure(xn,
                                                                                      previous_material.initial_density)

        # check tolerance condition
        if -tol < y1 < tol:
            new_intersection.isentrope = self.calculate_isentrope(hugoniot=previous_material_hugoniot,
                                                                  intersection=new_intersection,
                                                                  Gamma_eff=previous_material.Gamma_eff,
                                                                  initial_volume=previous_material.initial_volume,
                                                                  released=previous_material.released)
            new_intersection.isentrope.intersection=next_material_intersection
            return new_intersection, n

        # recursive call with updated interval
        return self.__secant_method(intersection1, new_intersection,
                                    previous_material_hugoniot,
                                    next_material_intersection,
                                    previous_material, next_material,
                                    n=n)
