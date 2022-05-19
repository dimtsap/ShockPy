from shock_wave_compression.Material import Material
from shock_wave_compression.Kapton import Kapton
from shock_wave_compression.MgO import MgO
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

from shock_wave_compression.Quartz import Quartz

# shock_pressure = 200  # GPa
# ############
# #  Step-1  #
# ############
# kapton = Kapton()
# up1, V1, compression_ratio = kapton.find_hugoniot_point_at_pressure(shock_pressure)
# (pressures_kapton, up_kapton, Us_kapton, volumes_kapton) = kapton.calculate_nominal_hugoniot()
# release_pressures, release_up = kapton.calculate_isentrope(shock_pressure, V1, up1, compression_ratio)
#
# mgo = MgO()
# (pressuresMgO, upMgO, UsMgO, volumesMgO) = mgo.calculate_nominal_hugoniot()
#
# # mgo_up, mgo_p = Material.calculate_intersection(upMgO, pressuresMgO, release_up[0], release_pressures[0])
# #
# # if up1 > mgo_up:
# #     indices = np.where((mgo_up <= release_up[0]) & (release_up[0] <= up1))
# # else:
# #     indices = np.where((up1 <= release_up[0]) & (release_up[0] <= mgo_up))
# # release_kapton_up = release_up[0][indices]
# # release_kapton_p = release_pressures[0][indices]
# #
# # ############
# # #  Step-2  #
# # ############
# # up2, V2, compression_ratio_2 = mgo.find_hugoniot_point_at_pressure(mgo_p)
# # (pressures_mgo, up_mgo, Us_mgo, volumes_mgo) = mgo.calculate_nominal_hugoniot()
# # release_pressures_MgO, release_up_MGO = mgo.calculate_isentrope(mgo_p, V2, up2, compression_ratio_2)
# #
# # quartz = Quartz()
# # (pressures_quartz, up_quartz, Us_quartz, volumes_quartz) = quartz.calculate_nominal_hugoniot()
# #
# # quartz_up, quartz_p = Material.calculate_intersection(up_quartz, pressures_quartz,
# #                                                       np.array(release_up_MGO[0][1:]),
# #                                                       np.array(release_pressures_MgO[0][1:]))
# #
# # if up2 > quartz_up:
# #     indices = np.where((quartz_up <= release_up_MGO[0]) & (release_up_MGO[0] <= up2))
# # else:
# #     indices = np.where((up2 <= release_up_MGO[0]) & (release_up_MGO[0] <= quartz_up))
# # release_up_MGO = release_up_MGO[0][indices]
# # release_pressures_MgO = release_pressures_MgO[0][indices]
# #
# # ###########
#
# plt.scatter(up1, shock_pressure, color='#af2c20')
# plt.plot(up_kapton, pressures_kapton, color="r", label='Kapton Hugoniot')
# # plt.plot(release_kapton_up, release_kapton_p, color="#e26c62", linestyle='dashed', label='Kapton isentrope')
# #
# # plt.plot(upMgO, pressuresMgO, color="b", label='MgO Hugoniot')
# # plt.plot(release_up_MGO, release_pressures_MgO, color="#ADD8E6", linestyle='dashed', label='MgO isentrope')
# # plt.scatter(mgo_up, mgo_p, color='#1c2e4a')
# #
# # plt.plot(up_quartz, pressures_quartz, color="black", label='Quartz Hugoniot')
# # plt.scatter(quartz_up, quartz_p, color="gray")
# plt.xlim((0, 12))
# plt.ylim((0, 700))
# plt.legend()
# plt.show()


from shock_wave_compression.ShockWaveExperiment import ShockWaveExperiment

shock_pressure = 200  # GPa

experiment = ShockWaveExperiment(materials=[Kapton(), MgO(), Quartz()])

experiment.run_experiment(shock_pressure)