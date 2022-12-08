import os.path

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

from shock_wave_compression.material_states.Hugoniot import Hugoniot
from shock_wave_compression.material_database.baseclass.Material import Material
from shock_wave_compression.material_states.Isentrope import Isentrope


class Kapton(Material):

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 1.33):  # Source https://dielectricmfg.com/knowledge-base/kapton/
        super().__init__(1.1, initial_density)
        hugoniot_P_up = np.loadtxt(os.path.join(os.getcwd(), "material_database", "data", "kapton_hugoniot.txt"))
        p = hugoniot_P_up[:, 0]
        up = hugoniot_P_up[:, 1]
        us = p / (self.initial_density * up)
        v = self.initial_volume * (us - up) / us
        self.nominal_hugoniot = Hugoniot(pressures=p, particle_velocities=up, shock_velocities=us, volumes=v)
        self.hugoniots_list = [self.nominal_hugoniot]


    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        pass
