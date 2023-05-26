import os.path

import numpy as np
from scipy import interpolate

from shock_wave_compression.material_states.Hugoniot import Hugoniot
from shock_wave_compression.material_database.baseclass.Material import Material
from scipy.integrate import quad, cumulative_trapezoid

from shock_wave_compression.material_states.Isentrope import Isentrope


class Kapton(Material):

    def __init__(self, gamma_eff: float = 0.7,
                 initial_density: float = 1.33):  # Source https://dielectricmfg.com/knowledge-base/kapton/
        super().__init__(gamma_eff, initial_density)
        hugoniot_P_up = np.loadtxt(os.path.join(os.getcwd(), "material_database", "data", "kapton_hugoniot.txt"))
        
        p = hugoniot_P_up[:, 0]
        up = hugoniot_P_up[:, 1]
        us = p / (self.initial_density * up)
        v = self.initial_volume * (us - up) / us
        self.nominal_hugoniot = Hugoniot(pressures=p, particle_velocities=up, shock_velocities=us, volumes=v)
        self.hugoniots_list = [self.nominal_hugoniot]

    def calculate_isentrope(self, hugoniot, intersection):
        pressuresList = hugoniot.pressures.tolist()
        volumesList = hugoniot.volumes.tolist()
        release_pressures = []
        release_particle_velocities = []
        V = np.linspace(intersection.volume, min(volumesList), num=1000)
        hugoniot_interpolator = interpolate.interp1d(np.flip(np.array(volumesList)), np.flip(np.array(pressuresList)),
                                                     kind='cubic',
                                                     fill_value="extrapolate")
        reference_hugoniot_pressure = hugoniot_interpolator(V)
        # for (volume_Hugoniot, pressure_Hugoniot) in zip([volumesList], [pressuresList]):
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

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        pass
