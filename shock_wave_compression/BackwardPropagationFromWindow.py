from shock_wave_compression.material_database.baseclass.Material import Material
import numpy as np
import seaborn as sns
import colorsys
import matplotlib.pyplot as plt


class BackwardPropagationFromWindow:

    def __init__(self, materials: list[Material]):
        self.initial_intersections = None
        self.materials = materials
        self.n_materials = len(materials)
        self._palette = sns.color_palette(None, self.n_materials)
        self._very_light_palette = [self._scale_lightness(color, 1.8) for color in self._palette]
        self._lighter_palette = [self._scale_lightness(color, 1.5) for color in self._palette]
        self._darker_palette = [self._scale_lightness(color, 0.5) for color in self._palette]

    @staticmethod
    def _scale_lightness(rgb, scale_l):
        # convert rgb to hls
        h, l, s = colorsys.rgb_to_hls(*rgb)
        # manipulate h, l, s values and return as rgb
        return colorsys.hls_to_rgb(h, min(1, l * scale_l), s=s)

    def propagate(self, measured_pressure: float, std=None):
        if std is not None:
            pressures_range = np.linspace((1 - 2 * std) * measured_pressure,
                                          (1 + 2 * std) * measured_pressure, num=5) \
                .tolist()
        else:
            pressures_range = [measured_pressure]

        intersections = []
        for pressure in pressures_range:
            intersections.extend(self.materials[-1].intersections_at_pressure(pressure))

        for index_material in range(self.n_materials - 2, -1, -1):
            material_i = self.materials[index_material]
            previous_material_intersections = []
            for intersection in intersections:
                previous_material_intersections.extend(material_i.find_hugoniots_from_next_intersection(intersection))

            intersections = previous_material_intersections

        self.initial_intersections = intersections

    def plot(self):
        fig, ax = plt.subplots(1, 1)
        for initial_intersection in self.initial_intersections:
            index_material = 0
            self._plot_intersection(initial_intersection, index_material, ax)
        for index_material in range(self.n_materials):
            material = self.materials[index_material]
            if material.is_stochastic:
                for index_hugoniot in range(len(material.hugoniots_list)):
                    hugoniot = material.hugoniots_list[index_hugoniot]
                    plt.plot(hugoniot.particle_velocities, hugoniot.pressures, color='gray', zorder=0)
            plt.plot(material.nominal_hugoniot.particle_velocities,
                     material.nominal_hugoniot.pressures, color=self._palette[index_material],
                     label=material.__class__.__name__ + " Hugoniot", zorder=1)

        plt.xlim((4, 12))
        plt.ylim((0, 700))
        plt.legend()
        plt.ylabel("Pressure (GPa)")
        plt.xlabel("Particle Velocity (km/s)")
        plt.savefig("uncertain_backward.png")

    def _plot_isentrope(self, isentrope, index_material, ax):
        ax.plot(np.squeeze(isentrope.particle_velocities),
                np.squeeze(isentrope.pressures),
                color=self._lighter_palette[index_material], linestyle='dashed', zorder=1)
        if isentrope.intersection is not None:
            index_material += 1
            self._plot_intersection(isentrope.intersection, index_material=index_material, ax=ax)

    def _plot_intersection(self, intersection, index_material, ax):
        ax.scatter(intersection.particle_velocity, intersection.pressure,
                   color=self._darker_palette[index_material], zorder=2)
        if intersection.isentrope is not None:
            self._plot_isentrope(intersection.isentrope, index_material, ax=ax)