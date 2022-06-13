import os.path

import numpy as np
from scipy import interpolate

from shock_wave_compression.material_states.Hugoniot import Hugoniot
from shock_wave_compression.material_database.baseclass.Material import Material


class Kapton(Material):

    def calculate_stochastic_hugoniot(self) -> list[Hugoniot]:
        pass

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 1.33):  # Source https://dielectricmfg.com/knowledge-base/kapton/
        super().__init__(1.1, initial_density)
        hugoniot_P_up = np.loadtxt(os.path.join(os.getcwd(), "material_database", "data", "kapton_hugoniot.txt"))
        # hugoniot_P_up = np.loadtxt('./data/kapton_hugoniot.txt')
        p = hugoniot_P_up[:, 0]
        up = hugoniot_P_up[:, 1]
        us = p / (self.initial_density * up)
        v = self.initial_volume * (us - up) / us
        self.nominal_hugoniot = Hugoniot(pressure=p, particle_velocity=up, shock_velocity=us, volume=v)
        self.hugoniots_list = [self.nominal_hugoniot]
        self.intersections = []

    def calculate_nominal_hugoniot(self):
        return self.nominal_hugoniot

    def calculate_isentrope(self, hugoniot, intersection):
        isentrope = super().calculate_isentrope(hugoniot, intersection)
        # tries to extrapolate Kapton isentrope
        f = interpolate.interp1d(isentrope.particle_velocities[0, 1:],
                                 isentrope.pressures[0, 1:],
                                 kind='slinear', fill_value='extrapolate')
        # release_up[0][1:], release_pressures[0][1:],
        release_kapton_up = np.linspace(2, 15)
        release_kapton_p = f(release_kapton_up)

        isentrope.pressures = release_kapton_p
        isentrope.particle_velocities = release_kapton_up
        return isentrope

    def release_isentropes_at_pressure(self, pressure: float):
        isentropes=[]
        for hugoniot in self.hugoniots_list:
            intersection = self._find_hugoniot_point_at_pressure(hugoniot, pressure)
            release_isentrope = self.calculate_isentrope(hugoniot, intersection)
            isentropes.append(release_isentrope)
        return isentropes
