from ShockPy.shock_wave_compression.material_states.Isentrope import Isentrope


class ReflectedHugoniot:

    def calculate_isentrope(self, hugoniot, intersection, Gamma_eff, initial_volume, released=True):
        """
        This method takes as input the Hugoniot and the current material shocked state as the intersection
        and mirrors at alog the :math:`u_p=u_{p_{shocked}}` to computed the reflected Hugoniot approximation of the
        isentrope.

        :param hugoniot: Hugoniot of the current material
        :param intersection: An :class:`.Intersection` object containing the information about the current material
         shocked state.
        :param Gamma_eff: This parameter is used to maintain a common interface between the isentrope calculators.
         Not required for the Reflected Hugoniot case.
        :param initial_volume: This parameter is used to maintain a common interface between the isentrope calculators.
         Not required for the Reflected Hugoniot case.
        :param released: This parameter is used to maintain a common interface between the isentrope calculators.
         Not required for the Reflected Hugoniot case.
        """
        reflected_hugoniot_particle_velocities = intersection.particle_velocity + \
                                                 (intersection.particle_velocity - hugoniot.particle_velocities)

        return Isentrope(
            pressures=hugoniot.pressures,
            particle_velocities=reflected_hugoniot_particle_velocities,
            intersection=intersection,
        )

    def find_previous_material_isentrope(self,
                                         previous_material_hugoniot,
                                         current_material_intersection,
                                         previous_material,
                                         next_material):
        """
        This methods evaluates the Isentrope of the previous material, that goes through the current material shocked
        state. The two material Hugoniots and the current material shocked state are required. This functions aids in
        the backward propagation of the experiment from the window to the input laser drive.

        :param previous_material_hugoniot: The Hugoniot of the previous material. In both Isentrope Calculators it is
         useful for obtaining the required isentrope and thus propagating the experiment backwards.
        :param current_material_intersection: An :class:`Intersection` object that contains the information of the
         current material shocked state.
        :param previous_material: An object reference to the previous material and its properties
        :param next_material: An object reference to the current material and its properties
        """
        previous_material_density = previous_material.initial_density
        negative_hugoniot = previous_material_hugoniot.reflected_hugoniot(previous_material_density)
        negative_hugoniot_intersection = \
            negative_hugoniot.find_hugoniot_point_at_pressure(current_material_intersection.pressure,
                                                              previous_material_density)
        delta_particle_velocity = current_material_intersection.particle_velocity - \
                                  negative_hugoniot_intersection.particle_velocity
        reflected_hugoniot = previous_material_hugoniot \
            .reflected_hugoniot(previous_material_density,
                                delta_particle_velocity=delta_particle_velocity)

        return Isentrope(
            pressures=reflected_hugoniot.pressures,
            particle_velocities=reflected_hugoniot.particle_velocities,
            intersection=current_material_intersection,
        )
