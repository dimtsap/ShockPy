from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
import numpy as np
import seaborn as sns
import colorsys
import matplotlib.pyplot as plt
# import scienceplots


class ShockWaveExperiment:

    def __init__(self, materials: list[Material]):
        """
        Class simulating the forward propagation of experiment using the same sequence of materials as in the
        experimental setup.

        :param materials: A list containing :class:`.Material` objects in the same sequence as in the experiment.
        """
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

    def run_experiment(self,
                       shock_pressure: float,
                       cov=None,
                       initial_points: int = 500,
                       n_isentropes_per_material: int = None):
        """
        Upon initialization of the class, this method is used to forward propagate the experiment using the Impedance
        Matching technique.

        :param shock_pressure: Predefined pressure of the input laser drive.
        :param cov: Coefficient of variation of pressure assumed for the laser drive. If value of this parameter is
         :any:`None`, the algorithm will assume there is no uncertainty in the input.
        :param initial_points: This parameter defines the number of points to be drawn from a Normal distribution
         with mean value the :code:`shock_pressure` and standard deviation inferred using the :code:`cov` parameter.
         This will create a number of initial shocked states where the experiment will start from, thus taking into
         account the uncertainty of the input laser drive.
        :param n_isentropes_per_material: This parameter limits the number of isentropes generated and propagated from
         one material to the next, and as a result prevents the number of experiment evaluations from growing
         geometrically.
        """
        pressures_range = [shock_pressure] if cov is None else \
            cov * shock_pressure * np.random.randn(initial_points) + shock_pressure

        release_isentropes = []
        for pressure in pressures_range:
            release_isentropes.extend(self.materials[0].release_isentropes_at_pressure(pressure))
        for index_material in range(1, self.n_materials):
            material_i = self.materials[index_material]

            new_material_isentropes = []
            for isentrope in release_isentropes:
                new_material_isentropes.extend(material_i.hugoniots_intersection_with_isentrope(isentrope))

            release_isentropes = new_material_isentropes if n_isentropes_per_material is None \
                else [new_material_isentropes[index] for index in
                      list(np.random.randint(0, len(new_material_isentropes),
                                             size=n_isentropes_per_material))]

        self.final_isentropes = release_isentropes

    def plot(self):
        """
        This functions plots all possible experimental paths the forward propagation of an experiment has taken into
        account, creating a plot similar to the one below.
        """
        # plt.style.use('science')

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

        from datetime import datetime

        plt.xlim((0, 12))
        plt.ylim((0, 700))
        plt.ylabel("Pressure (GPa)", fontsize=18)
        plt.xlabel("Particle Velocity (km/s)", fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        hyades_up = [9.65, 6.81, 8.15]
        hyades_p = [211.39, 373.47, 304.85]
        # plt.axhline(y=211.39, color='r', linestyle='-')
        # plt.axhline(y=373.47, color='r', linestyle='-')
        # plt.axhline(y=304.85, color='r', linestyle='-')
        # plt.scatter(x=hyades_up[0], y=hyades_p[0], color="black", marker='*', label='HYADES Kapton point', s=160,
        #             zorder=1e6)
        # plt.scatter(x=hyades_up[1], y=hyades_p[1], color="black", marker='*', label='HYADES MgO point', s=160,
        #             zorder=1e6)
        # plt.scatter(x=hyades_up[2], y=hyades_p[2], color="black", marker='*', label='HYADES Quartz point', s=160,
        #             zorder=1e6)
        # plt.scatter(x=hyades_up, y=hyades_p, c='black', marker='*')
        # plt.scatter(x=[x + 9.33 - 4.54 for x in hyades_up], y=hyades_p, c='black', marker='x')
        # plt.legend(prop={'size': 16}, loc='upper left')
        # plt.savefig(f"uncertain_shockwave_{datetime.now()}.png")
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
