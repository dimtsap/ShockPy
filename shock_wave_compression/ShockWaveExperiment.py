from shock_wave_compression.material_database.baseclass.Material import Material
import numpy as np
import seaborn as sns
import colorsys
import matplotlib.pyplot as plt


class ShockWaveExperiment:

    def __init__(self, materials: list[Material]):
        self.final_isentropes = None
        self.materials = materials
        self.n_materials = len(self.materials)
        self._palette = sns.color_palette(None, self.n_materials)
        self._very_light_palette = [ShockWaveExperiment._scale_lightness(color, 1.8) for color in self._palette]
        self._lighter_palette = [ShockWaveExperiment._scale_lightness(color, 1.5) for color in self._palette]
        self._darker_palette = [ShockWaveExperiment._scale_lightness(color, 0.5) for color in self._palette]

    @staticmethod
    def _scale_lightness(rgb, scale_l):
        # convert rgb to hls
        h, l, s = colorsys.rgb_to_hls(*rgb)
        # manipulate h, l, s values and return as rgb
        return colorsys.hls_to_rgb(h, min(1, l * scale_l), s=s)

    def run_experiment(self, shock_pressure: float, std=None):
        if std is not None:
            pressures_range = np.linspace((1-2*std)*shock_pressure, (1+2*std)*shock_pressure, num=5).tolist()
        else:
            pressures_range = [shock_pressure]
        release_isentropes=[]
        for pressure in pressures_range:
            release_isentropes.extend(self.materials[0].release_isentropes_at_pressure(pressure))
        # release_isentropes =
        for index_material in range(1, self.n_materials):
            material_i = self.materials[index_material]

            new_material_isentropes = []
            for isentrope in release_isentropes:
                new_material_isentropes.extend(material_i.hugoniots_intersection_with_isentrope(isentrope))

            release_isentropes = new_material_isentropes

        self.final_isentropes = release_isentropes

    def plot(self):
        plt.style.use('science')

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        for final_isentrope in self.final_isentropes:
            index_material = len(self.materials) - 1
            self._plot_intersection(final_isentrope.intersection, index_material, ax)
            # self._plot_isentrope(isentrope=final_isentrope, index_material=index_material, ax=ax)
        for index_material in range(self.n_materials):
            material = self.materials[index_material]
            if material.is_stochastic:
                for index_hugoniot in range(len(material.hugoniots_list)):
                    hugoniot = material.hugoniots_list[index_hugoniot]
                    plt.plot(hugoniot.particle_velocities, hugoniot.pressures, color='gray', linewidth=2, zorder=0)
            plt.plot(material.nominal_hugoniot.particle_velocities,
                     material.nominal_hugoniot.pressures, color=self._palette[index_material],
                     label=material.__class__.__name__ + " Hugoniot", linewidth=2, zorder=1)

        plt.xlim((4, 12))
        plt.ylim((0, 700))
        plt.legend(prop={'size': 16}, loc='upper left')
        plt.ylabel("Pressure (GPa)", fontsize=18)
        plt.xlabel("Particle Velocity (km/s)", fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig("uncertain_shockwave.png")
        plt.show()

    def _plot_isentrope(self, isentrope, index_material, ax):
        ax.plot(np.squeeze(isentrope.particle_velocities),
                np.squeeze(isentrope.pressures),
                  color=self._lighter_palette[index_material], linestyle='dashed', zorder=1)
        if isentrope.intersection is not None:
            self._plot_intersection(isentrope.intersection, index_material=index_material, ax=ax)

    def _plot_intersection(self, intersection, index_material, ax):
        ax.scatter(intersection.particle_velocity, intersection.pressure,
                     color=self._darker_palette[index_material], zorder=2)
        if intersection.isentrope is not None:
            index_material -= 1
            self._plot_isentrope(intersection.isentrope, index_material, ax=ax)
