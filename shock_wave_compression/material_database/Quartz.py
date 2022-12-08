import math

import numpy as np
from scipy import interpolate
from shock_wave_compression.material_database.baseclass.Material import Material
from shock_wave_compression.material_states.Hugoniot import Hugoniot
import os
import pickle

from shock_wave_compression.material_states.Isentrope import Isentrope


class Quartz(Material):

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 2.65,
                 is_stochastic=False):
        super().__init__(gamma_eff, initial_density, is_stochastic)
        with open(os.path.join(os.getcwd(), "material_database", "data", "samples_quartz_exponential0.p"), "rb") \
                as input_file:
            data = pickle.load(input_file)
        # self.nominal_hugoniot = self.calculate_hugoniot(np.array([1.754, 1.862, -3.364e-2, 5.666e-4])) # polynomial hugoniot
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.278, 1.193, -2.505, -0.3701]))  # exponential hugoniot
        if is_stochastic:
            self.hugoniots_list = self.hugoniots_list = [self.calculate_hugoniot(x) for x in data[:1000:100]]
        else:
            self.hugoniots_list = [self.nominal_hugoniot]

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        return np.array([parameters[0] + parameters[1] * x + parameters[2] * x * math.exp(parameters[3] * x)
                         for x in hugoniot_particle_velocity])
