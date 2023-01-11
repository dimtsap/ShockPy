from shock_wave_compression.material_database.Kapton import Kapton
from shock_wave_compression.material_database.MgO import MgO
from shock_wave_compression.material_database.Quartz import Quartz


from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

shock_pressure = 200  # GPa

experiment = ShockWaveExperiment(materials=[Kapton(), MgO(is_stochastic=False), Quartz(is_stochastic=False)])

experiment.run_experiment(shock_pressure)

experiment.plot()

