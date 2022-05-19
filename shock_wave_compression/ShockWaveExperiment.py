from shock_wave_compression.Material import Material
import numpy as np
import seaborn as sns
import colorsys
import matplotlib.pyplot as plt


class ShockWaveExperiment:

    def __init__(self, materials: list[Material]):
        self.materials = materials
        self.n_materials = len(self.materials)
        self._palette = sns.color_palette(None, self.n_materials)
        self._lighter_palette = [ShockWaveExperiment._scale_lightness(color, 1.5) for color in self._palette]
        self._darker_palette = [ShockWaveExperiment._scale_lightness(color, 0.5) for color in self._palette]

    @staticmethod
    def _scale_lightness(rgb, scale_l):
        # convert rgb to hls
        h, l, s = colorsys.rgb_to_hls(*rgb)
        # manipulate h, l, s values and return as rgb
        return colorsys.hls_to_rgb(h, min(1, l * scale_l), s=s)

    def run_experiment(self, shock_pressure: float):
        input_pressure = shock_pressure
        for index_material in range(self.n_materials - 1):
            material_i = self.materials[index_material]

            up1, V1, compression_ratio = material_i.find_hugoniot_point_at_pressure(input_pressure)
            plt.scatter(up1, input_pressure, color=self._darker_palette[index_material])

            (pressures_hugoniot_i, up_hugoniot_i,
             Us_hugoniot_i, volumes_hugoniot_i) = material_i.calculate_nominal_hugoniot()
            plt.plot(up_hugoniot_i, pressures_hugoniot_i, color=self._palette[index_material],
                     label=material_i.__class__.__name__ + " Hugoniot")

            release_p, release_up = material_i.calculate_isentrope(input_pressure, V1, up1, compression_ratio)

            material_i_1 = self.materials[index_material+1]
            (pressures_i_1, up_i_1, Us_i_1, volumes_i_1) = material_i_1.calculate_nominal_hugoniot()

            intersection_up, intersection_p = \
                Material.calculate_intersection(up_i_1, pressures_i_1, release_up[0], release_p[0])
            plt.scatter(intersection_up, intersection_p, color=self._darker_palette[index_material+1])

            if up1 > intersection_up:
                indices = np.where((intersection_up <= release_up[0]) & (release_up[0] <= up1))
            else:
                indices = np.where((up1 <= release_up[0]) & (release_up[0] <= intersection_up))
            release_up = release_up[0][indices]
            release_p = release_p[0][indices]
            plt.plot(release_up, release_p, color=self._lighter_palette[index_material], linestyle='dashed')

            input_pressure = intersection_p

        material_last = self.materials[-1]
        (pressures_last, up_last, Us_last, volumes_last) = material_last.calculate_nominal_hugoniot()
        plt.plot(up_last, pressures_last, color=self._palette[-1],
                 label=material_last.__class__.__name__ + " Hugoniot")

        plt.xlim((4, 12))
        plt.ylim((0, 700))
        plt.legend()
        plt.ylabel("Pressure (GPa)")
        plt.xlabel("Particle Velocity (km/s)")
        plt.show()
