import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import interpolate

from shock_wave_compression.material_states.Isentrope import Isentrope


class IntegratedIsentrope:

    def calculate_isentrope(self, hugoniot, intersection, Gamma_eff, initial_volume, released=True):
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