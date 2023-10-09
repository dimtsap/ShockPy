"""

Uncertain forward shock-wave experiment propagation
=========================================================
"""

from ShockPy.shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment
from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
from ShockPy.shock_wave_compression.material_states.isentrope_calculators.ReflectedHugoniot import \
    ReflectedHugoniot

from ShockPy.shock_wave_compression.material_database.Kapton import Kapton
from ShockPy.shock_wave_compression.material_database.MgO import MgO
from ShockPy.shock_wave_compression.material_database.Quartz import Quartz

shock_pressure = 211.39  # GPa

Material.isentrope_calculator = ReflectedHugoniot()

experiment = ShockWaveExperiment(materials=[Kapton(released=False),
                                            MgO(is_stochastic=True, released=True),
                                            Quartz(is_stochastic=True, released=True)])
experiment.run_experiment(shock_pressure, cov=0.05, initial_points=100)

experiment.plot()

