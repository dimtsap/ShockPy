import os.path

import numpy as np
from src.shock_wave_compression.material_database.baseclass.Material import Material
import pickle


class MgO(Material):

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 3.58,
                 is_stochastic=False,
                 released=True):
        super().__init__(gamma_eff=gamma_eff,
                         initial_density=initial_density,
                         is_stochastic=is_stochastic,
                         released=released)
        with open(os.path.join(os.path.dirname(__file__), "data", "samplesMgO.p"), "rb") as input_file:
            data = pickle.load(input_file)
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.6161, 1.4111, -0.016277]))
        self.hugoniots_list = [self.nominal_hugoniot] if not is_stochastic \
            else [self.calculate_hugoniot(x) for x in data[:1000:10]]

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        return parameters[0] + parameters[1] * hugoniot_particle_velocity + parameters[2] * np.square(
            hugoniot_particle_velocity)