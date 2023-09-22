from ImpedancePy.shock_wave_compression.material_states.Isentrope import Isentrope


class ReflectedHugoniot:

    def calculate_isentrope(self, hugoniot, intersection, Gamma_eff, initial_volume, released=True):
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
