import math

import numpy as np
from src.shock_wave_compression.material_database.baseclass.Material import Material
import os
import pickle


class Quartz(Material):

    def __init__(self, gamma_eff: float = 0.6,
                 initial_density: float = 2.65,
                 released=True,
                 is_stochastic=False):
        super().__init__(gamma_eff, initial_density, is_stochastic, released)
        with open(os.path.join(os.path.dirname(__file__), "data", "samples_quartz_exponential0.p"), "rb") \
                as input_file:
            data = pickle.load(input_file)
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.278, 1.193, -2.505, -0.3701]))  # exponential hugoniot
        self.hugoniots_list = [self.nominal_hugoniot] if not is_stochastic \
            else [self.calculate_hugoniot(x) for x in data[:1000:10]]

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        return np.array([parameters[0] + parameters[1] * x + parameters[2] * x * math.exp(parameters[3] * x)
                         for x in hugoniot_particle_velocity])

