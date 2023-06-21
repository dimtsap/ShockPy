from shock_wave_compression.material_states.Isentrope import Isentrope


class ReflectedHugoniot:

    def calculate_isentrope(self, hugoniot, intersection, Gamma_eff, initial_volume, released=True):
        reflected_hugoniot_particle_velocities = intersection.particle_velocity + \
                                                 (intersection.particle_velocity - hugoniot.particle_velocities)

        return Isentrope(
            pressures=hugoniot.pressures,
            particle_velocities=reflected_hugoniot_particle_velocities,
            intersection=intersection,
        )
