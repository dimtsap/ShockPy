from shock_wave_compression.material_database.Kapton import Kapton
from shock_wave_compression.material_database.MgO import MgO
from shock_wave_compression.material_database.Quartz import Quartz
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

shock_pressure = 211.39  # GPa

experiment = ShockWaveExperiment(materials=[Kapton(), MgO(is_stochastic=False), Quartz(is_stochastic=False)])

experiment.run_experiment(shock_pressure)

experiment.plot()


window_intersections = np.array(
    [(isentrope.intersection.particle_velocity, isentrope.intersection.pressure) for isentrope in
     experiment.final_isentropes])
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