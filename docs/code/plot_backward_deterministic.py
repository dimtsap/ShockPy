"""

Deterministic backward shock-wave experiment propagation
=========================================================
"""

from ShockPy.shock_wave_compression.material_database.baseclass import Material
from ShockPy.shock_wave_compression.material_states.isentrope_calculators import ReflectedHugoniot, \
    IntegratedIsentrope
from ShockPy.shock_wave_compression.BackwardPropagationFromWindow import BackwardPropagationFromWindow
from ShockPy.shock_wave_compression.material_database.Kapton import Kapton
from ShockPy.shock_wave_compression.material_database.MgO import MgO
from ShockPy.shock_wave_compression.material_database.Quartz import Quartz

Material.isentrope_calculator = IntegratedIsentrope()

measured_pressure = 315.03  # GPa
backward = BackwardPropagationFromWindow(materials=[Kapton(released=False),
                                                    MgO(is_stochastic=False, released=True),
                                                    Quartz(is_stochastic=False, released=True)])
backward.propagate(measured_pressure)
backward.plot()

