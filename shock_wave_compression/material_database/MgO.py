import os.path

import numpy as np
from scipy import interpolate
from shock_wave_compression.material_database.baseclass.Material import Material
from shock_wave_compression.material_states.Hugoniot import Hugoniot
import pickle


class MgO(Material):

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 3.58,
                 is_stochastic=False):
        super().__init__(gamma_eff, initial_density, is_stochastic)
        with open(os.path.join(os.getcwd(), "material_database", "data", "samplesMgO.p"), "rb") as input_file:
            data = pickle.load(input_file)
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.6161, 1.4111, -0.016277]))
        if is_stochastic:
            self.hugoniots_list = [self.calculate_hugoniot(x) for x in data[:1000:20]]
        else:
            self.hugoniots_list = [self.nominal_hugoniot]

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        return parameters[0] + parameters[1] * hugoniot_particle_velocity +parameters[2] * np.square(
            hugoniot_particle_velocity)