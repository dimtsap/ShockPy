from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
import numpy as np
import seaborn as sns
import colorsys
import matplotlib.pyplot as plt
# import scienceplots


class BackwardPropagationFromWindow:

    def __init__(self, materials: list[Material]):
        """
        Class simulating the backward propagation of experiment using the same sequence of materials as in the
        experimental setup.

        :param materials: A list containing :class:`.Material` objects in the same sequence as in the experiment.
        """
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

    def propagate(self,
                  measured_pressure: float,
                  cov=None,
                  initial_points=1000,
                  n_intersections_per_material: int = None):
        """
        Upon initialization of the class, this method is used to backward tracethe experiment using the Impedance
        Matching technique.

        :param measured_pressure: Pressure measure at the window material layer of the experimental setup.
        :param cov: Coefficient of variation of pressure measured at the window layer. If value of this parameter is
         :any:`None`, the algorithm will assume there is no uncertainty in the measurement.
        :param initial_points: This parameter defines the number of points to be drawn from a Normal distribution
         with mean value the :code:`measured_pressure` and standard deviation inferred using the :code:`cov` parameter.
         This will create a number of initial measured pressures where the experiment will start from, thus taking into
         account the uncertainty of the window measurement.
        :param n_intersections_per_material: This parameter limits the number of shocked states calculated during the
         backward propagation of the experiment, and as a result prevents the number of experiment evaluations from
         growing geometrically.
        """

        pressures_range = [measured_pressure] if cov is None else \
            cov * measured_pressure * np.random.randn(initial_points) + measured_pressure

        intersections = []
        for index_pressure, pressure in enumerate(pressures_range):
            print(f'Material: {self.materials[-1].__class__}, pressure: {index_pressure}')
            intersections.extend(self.materials[-1].intersections_at_pressure(pressure))

        for index_material in range(self.n_materials - 2, -1, -1):
            previous_material_intersections = intersections if n_intersections_per_material is None \
                else [intersections[index] for index in
                      list(np.random.randint(0, len(intersections), size=n_intersections_per_material))]

            material_next = self.materials[index_material+1]
            material_i = self.materials[index_material]
            current_material_intersections = []
            for index_intersection, intersection in enumerate(previous_material_intersections):
                print(f'Material: {material_i.__class__}, intersection: {index_intersection}')
                current_material_intersections.extend(material_i
                                                      .find_previous_material_intersections_from_current(intersection,
                                                                                                         material_next))
            intersections = current_material_intersections

        self.initial_intersections = intersections

    def plot(self):
        """
        This functions plots all possible experimental paths the backward propagation of an experiment has taken into
        account, creating a plot similar to the one below.
        """
        # plt.style.use('science')
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
