"""

Deterministic forward shock-wave experiment propagation
=========================================================
"""

from src.ImpedancePy.shock_wave_compression.material_database.Kapton import Kapton
from src.ImpedancePy.shock_wave_compression.material_database.MgO import MgO
from src.ImpedancePy.shock_wave_compression.material_database.Quartz import Quartz
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from src.ImpedancePy.shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

# quartz = Quartz(is_stochastic=False)
# intersection = quartz._find_hugoniot_point_at_shock_velocity(quartz.nominal_hugoniot, 24.50)
from src.ImpedancePy.shock_wave_compression.material_database.baseclass.Material import Material
from src.ImpedancePy.shock_wave_compression.material_states.isentrope_calculators.ReflectedHugoniot import ReflectedHugoniot

shock_pressure = 211.39  # GPa

Material.isentrope_calculator = ReflectedHugoniot()

experiment = ShockWaveExperiment(materials=[Kapton(released=False),
                                            MgO(is_stochastic=False, released=True),
                                            Quartz(is_stochastic=False, released=True)])
experiment.run_experiment(shock_pressure)

experiment.plot()

