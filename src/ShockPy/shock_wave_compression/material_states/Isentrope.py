from dataclasses import dataclass
import numpy as np


@dataclass
class Isentrope:
    """
    Dataclass containing all info (i.e. :math:`P, u_p`) of an Isentrope.

    :param pressures: A numpy array containing the pressures :math:`(P)` of an Isentrope.
    :param particle_velocities: A numpy array containing the particle velocities :math:`(u_p)` of an Isentrope.
    :param intersection: An :class:`.Intersection` object define the shocked state of the material the initiated the
     isentrope.
    """
    pressures: np.ndarray
    particle_velocities: np.ndarray
    intersection: object
