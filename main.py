from shock_wave_compression.material_database.Kapton import Kapton
from shock_wave_compression.material_database.MgO import MgO
from shock_wave_compression.material_database.Quartz import Quartz
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

# quartz = Quartz(is_stochastic=False)
# intersection = quartz._find_hugoniot_point_at_shock_velocity(quartz.nominal_hugoniot, 24.50)
from shock_wave_compression.material_database.baseclass.Material import Material
from shock_wave_compression.material_states.isentrope_calculators.ReflectedHugoniot import ReflectedHugoniot

shock_pressure = 211.39  # GPa
# shock_pressure = intersection.pressure
# experiment = ShockWaveExperiment(materials=[quartz, MgO(is_stochastic=False)])
# Material.isentrope_calculator=ReflectedHugoniot()
experiment=ShockWaveExperiment(materials=[Kapton(released=False),
                                          MgO(is_stochastic=False,
                                              released=True),
                                          Quartz(is_stochastic=False,
                                                 released=True)])
experiment.run_experiment(shock_pressure, cov=0.05, initial_points=100)

experiment.plot()

window_points = [isentrope.intersection for isentrope in experiment.final_isentropes]
mgo_isentrope = [point.isentrope for point in window_points]
mgo_intersections = np.array(
    [(isentrope.intersection.particle_velocity, isentrope.intersection.pressure) for isentrope in
     mgo_isentrope])

window_intersections = np.array(
    [(isentrope.intersection.particle_velocity, isentrope.intersection.pressure) for isentrope in
     experiment.final_isentropes])
# gamma_eff = quartz.Gamma_eff
mean_window_pressure = np.mean(window_intersections[:, 1], dtype=np.float64)
std_window_pressure = np.std(window_intersections[:, 1], dtype=np.float64)
print(f"Mean window pressure:{mean_window_pressure}")
print(f"Std window pressure:{std_window_pressure}")

plt.style.use('science')
import datetime

df = pd.DataFrame(window_intersections, columns=['Particle Velocity (km/s)', 'Pressure (GPa)'])
sns.displot(data=df, x="Particle Velocity (km/s)", kind="kde")
plt.savefig(f'random_up_kde_{datetime.datetime.now()}.png')
sns.displot(data=df, x="Pressure (GPa)", kind="kde")
plt.savefig(f'random_P_kde_{datetime.datetime.now()}.png')
sns.displot(df, x="Particle Velocity (km/s)", y="Pressure (GPa)", kind='kde')
hyades_up = [9.65, 6.81, 8.15]
hyades_p = [211.39, 373.47, 304.85]
plt.scatter(x=hyades_up[2], y=hyades_p[2], color='black', marker='*', label='HYADES Quartz point', s=160)
plt.savefig(f'random_2D_kde_{datetime.datetime.now()}.png')

a = 1
