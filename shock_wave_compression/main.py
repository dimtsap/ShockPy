from shock_wave_compression.BackwardPropagationFromWindow import BackwardPropagationFromWindow
from shock_wave_compression.material_database.Kapton import Kapton
from shock_wave_compression.material_database.MgO import MgO
from shock_wave_compression.material_database.Quartz import Quartz


from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

# shock_pressure = 200  # GPa
#
# experiment = ShockWaveExperiment(materials=[Kapton(), MgO(is_stochastic=False), Quartz(is_stochastic=False)])
#
# experiment.run_experiment(shock_pressure)
#
# experiment.plot()

measured_pressure = 200  # GPa

backward = BackwardPropagationFromWindow(materials=[Kapton(), MgO(is_stochastic=False), Quartz(is_stochastic=False)])

backward.propagate(measured_pressure)
backward.plot()

