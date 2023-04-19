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

    def propagate(self, measured_pressure: float, std=None, initial_points=1000):
        if std is not None:
            # pressures_range = 2*1.96*std*measured_pressure*np.random.rand(initial_points)+(measured_pressure-1.96*std*measured_pressure)
            pressures_range = std * measured_pressure * np.random.randn(initial_points) + measured_pressure
            # pressures_range = np.linspace((1 - 1.96 * std) * measured_pressure, (1 + 1.96 * std) * measured_pressure,
            #                               num=initial_points).tolist()
        else:
            pressures_range = [measured_pressure]

        intersections = []
        for index_pressure, pressure in enumerate(pressures_range):
            print(f'Material: Quartz, pressure: {index_pressure}')
            intersections.extend(self.materials[-1].intersections_at_pressure(pressure))

        for index_material in range(self.n_materials - 2, -1, -1):
            previous_material_intersections = intersections
            intersections_1K = list(np.random.randint(0, len(previous_material_intersections), size=1000))
            previous_material_intersections = [previous_material_intersections[index] for index in intersections_1K]
            material_i = self.materials[index_material]
            current_material_intersections = []
            for index_intersection, intersection in enumerate(previous_material_intersections):
                print(f'Material: {material_i.__class__}, intersection: {index_intersection}')
                current_material_intersections.extend(material_i.find_hugoniots_from_next_intersection(intersection))
            previous_material_intersections = current_material_intersections

        self.initial_intersections = previous_material_intersections

    def plot(self):
        plt.style.use('science')
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        for initial_intersection in self.initial_intersections:
            index_material = 0
            self._plot_intersection(initial_intersection, index_material, ax)
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
        plt.savefig("uncertain_backward.png")
        plt.show()

    def _plot_isentrope(self, isentrope, index_material, ax):
        ax.plot(np.squeeze(isentrope.particle_velocities),
                np.squeeze(isentrope.pressures),
                color=self._lighter_palette[index_material], linewidth=2, linestyle='dashed', zorder=1)
        if isentrope.intersection is not None:
            index_material += 1
            self._plot_intersection(isentrope.intersection, index_material=index_material, ax=ax)

    def _plot_intersection(self, intersection, index_material, ax):
        ax.scatter(intersection.particle_velocity, intersection.pressure,
                   color=self._darker_palette[index_material], s=35, zorder=2)
        if intersection.isentrope is not None:
            self._plot_isentrope(intersection.isentrope, index_material, ax=ax)