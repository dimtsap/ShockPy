from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import splprep, splev
from scipy.spatial import cKDTree
from scipy import interpolate
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from shapely.geometry import LineString
from scipy.interpolate import interpolate


class Material(ABC):

    def __init__(self, gamma_eff: float,
                 initial_density: float):
        self.initial_density = initial_density
        self.initial_volume = 1 / initial_density
        self.Gamma_eff = gamma_eff

    @abstractmethod
    def find_hugoniot_point_at_pressure(self, pressure: float):
        """
        If the hugoniot is deterministic should return only its values.
        Otherwise, it should return the minimum and maximum particle velocities at this pressure.

        :param pressure:
        """
        pass

    @abstractmethod
    def calculate_nominal_hugoniot(self):
        pass

    def calculate_isentrope(self, input_pressure,
                            input_volume,
                            input_particle_velocity,
                            compression_ratio):
        (pressuresList, upList, UsList, volumesList) = self.calculate_nominal_hugoniot()
        release_pressures = []
        release_particle_velocities = []
        for (volume_Hugoniot, pressure_Hugoniot) in zip([volumesList], [pressuresList]):
            V = np.linspace(input_volume, max(volume_Hugoniot), num=1000)
            # V = np.linspace(input_volume, 0.1, num=1000)

            S = 1.31867  # What is this parameter.? Is it material dependent.?
            C01 = input_pressure / (self.initial_density * input_particle_velocity) - S * input_particle_velocity

            reference_hugoniot_pressure = self.initial_density * C01 ** 2 * (self.initial_volume / V - 1) \
                                          * (self.initial_volume / V) / np.square(S - (S - 1) * self.initial_volume / V)

            energy_integral = []
            for i in range(len(V)):
                delta_energy_integral_function = lambda x: (x / input_volume) ** self.Gamma_eff * \
                                                           reference_hugoniot_pressure[i] \
                                                           * (1 - self.Gamma_eff / 2 * (self.initial_volume / x - 1))
                aux = quad(delta_energy_integral_function, input_volume, V[i])
                energy_integral.append(aux[0])

            Es_E0 = [input_pressure * self.initial_volume / 2 * (compression_ratio - 1) / compression_ratio * (
                    input_volume / x) ** self.Gamma_eff - (input_volume / x) ** self.Gamma_eff * y
                     for (x, y) in zip(V, energy_integral)]

            release_pressure = [PH * (1 - self.Gamma_eff / 2 * (self.initial_volume / v - 1)) + self.Gamma_eff / v * dE
                                for (PH, v, dE) in zip(reference_hugoniot_pressure, V, Es_E0)]

            release_particle_velocity = [float(input_particle_velocity)]

            ups = list(np.real(
                input_particle_velocity + cumulative_trapezoid(np.sqrt(-np.diff(release_pressure) / np.diff(V)), V[1:],
                                                               initial=0)))
            release_particle_velocity.extend(ups)

            release_pressures.append(np.array(release_pressure))
            release_particle_velocities.append(np.array(release_particle_velocity))

        return release_pressures, release_particle_velocities

    @staticmethod
    def calculate_intersection(up1, p1, up2, p2):
        line_1 = LineString(np.column_stack((up1, p1)))
        line_2 = LineString(np.column_stack((up2, p2)))
        intersection = line_1.intersection(line_2)
        return intersection.xy[0][0], intersection.xy[1][0]


    @staticmethod
    def _upsample_coords(coord_list):
        # s is smoothness, set to zero
        # k is degree of the spline. setting to 1 for linear spline
        tck, u = splprep(coord_list, k=1, s=0.0)
        upsampled_coords = splev(np.linspace(0, 1, 100), tck)
        return upsampled_coords
