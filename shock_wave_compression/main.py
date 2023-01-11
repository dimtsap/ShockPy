from shock_wave_compression.material_database.Kapton import Kapton
from shock_wave_compression.material_database.MgO import MgO
from shock_wave_compression.material_database.Quartz import Quartz
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

shock_pressure = 200  # GPa

experiment = ShockWaveExperiment(materials=[Kapton(), MgO(is_stochastic=False), Quartz(is_stochastic=False)])

experiment.run_experiment(shock_pressure)

experiment.plot()


window_intersections = np.array(
    [(isentrope.intersection.particle_velocity, isentrope.intersection.pressure) for isentrope in
     experiment.final_isentropes])
mean_window_pressure = np.mean(window_intersections[:, 1], dtype=np.float64)
std_window_pressure = np.std(window_intersections[:, 1], dtype=np.float64)

plt.style.use('science')
df = pd.DataFrame(window_intersections, columns=['Particle Velocity (km/s)', 'Pressure (GPa)'])
sns.displot(data=df, x="Particle Velocity (km/s)", kind="kde")
sns.displot(data=df, x="Pressure (GPa)", kind="kde")
sns.displot(df, x="Particle Velocity (km/s)", y="Pressure (GPa)", kind='kde')
a=1