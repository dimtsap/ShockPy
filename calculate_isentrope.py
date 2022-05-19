# MgO
import numpy as np
from scipy import interpolate
from scipy.integrate import quad, cumulative_trapezoid

wm = [6.6161, 1.4111, -0.016277]
initial_density_MgO = 3.584
particle_velocity_MgO = list(np.linspace(3.5, 10, 131))
shock_velocity_MgO = [wm[0] + wm[1] * x + wm[2] * x ** 2 for x in particle_velocity_MgO]
pressure_MgO = [initial_density_MgO * x * y for (x, y) in zip(shock_velocity_MgO, particle_velocity_MgO)]
initial_volume_MgO = 1 / initial_density_MgO
volume_MgO = [initial_volume_MgO * (x - y) / x for (x, y) in zip(shock_velocity_MgO, particle_velocity_MgO)]

Gamma_eff = 1

#####
particle_velocity_Hugoniot = list(np.linspace(1, 10, 100))
shock_velocity_Hugoniot = [wm[0] + wm[1] * x + wm[2] * x ** 2 for x in particle_velocity_Hugoniot]
volume_Hugoniot = [initial_volume_MgO * (x - y) / x for (x, y) in zip(shock_velocity_Hugoniot,
                                                                      particle_velocity_Hugoniot)]
pressure_Hugoniot = [initial_density_MgO * x * y for (x, y) in zip(shock_velocity_Hugoniot, particle_velocity_Hugoniot)]

#####

####### The following should be input variables
input_pressure = pressure_MgO[0]
input_volume = volume_MgO[0]
input_particle_velocity = particle_velocity_MgO[0]
# input_shock_velocity = shock_velocity_MgO[0]
compression_ratio = initial_volume_MgO / input_volume
###########

V = list(np.linspace(input_volume, 0.244, 1000))
interpolator = interpolate.interp1d(volume_Hugoniot, pressure_Hugoniot, kind='cubic')

reference_hugoniot_pressure = interpolator(V)
energy_integral = []
for i in range(len(V)):
    delta_energy_integral_function = lambda x: (x / input_volume) ** Gamma_eff * reference_hugoniot_pressure[i] \
                                               * (1 - Gamma_eff / 2 * (initial_volume_MgO / x - 1))
    aux = quad(delta_energy_integral_function, input_volume, V[i])
    energy_integral.append(aux[0])

Es_E0 = [input_pressure * initial_volume_MgO / 2 * (compression_ratio - 1) / compression_ratio * (
            input_volume/x)**Gamma_eff-(input_volume/x)**Gamma_eff*y
         for (x, y) in zip(V, energy_integral)]

release_pressure = [PH * (1 - Gamma_eff / 2 * (initial_volume_MgO / v - 1)) + Gamma_eff / v * dE
                    for (PH, v, dE) in zip(reference_hugoniot_pressure, V, Es_E0)]

release_particle_velocity = [input_particle_velocity]

ups = list(np.real(input_particle_velocity + cumulative_trapezoid(np.sqrt(-np.diff(release_pressure)/np.diff(V)), V[1:], initial=0)))
release_particle_velocity.extend(ups)
