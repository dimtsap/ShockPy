import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import cKDTree
from scipy import interpolate

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.axis('off')

def upsample_coords(coord_list):
    # s is smoothness, set to zero
    # k is degree of the spline. setting to 1 for linear spline
    tck, u = interpolate.splprep(coord_list, k=1, s=0.0)
    upsampled_coords = interpolate.splev(np.linspace(0, 1, 100), tck)
    return upsampled_coords

# target line
x_targ = [1, 2, 3, 4, 5, 6, 7, 8]
y_targ = [20, 100, 50, 120, 55, 240, 50, 25]
z_targ = [20, 100, 50, 120, 55, 240, 50, 25]
targ_upsampled = upsample_coords([x_targ, y_targ, z_targ])
targ_coords = np.column_stack(targ_upsampled)

# KD-tree for nearest neighbor search
targ_kdtree = cKDTree(targ_coords)

# line two
x2 = [3,4,5,6,7,8,9]
y2 = [25,35,14,67,88,44,120]
z2 = [25,35,14,67,88,44,120]
l2_upsampled = upsample_coords([x2, y2, z2])
l2_coords = np.column_stack(l2_upsampled)

# plot both lines
ax.plot(x_targ, y_targ, z_targ, color='black', linewidth=0.5)
ax.plot(x2, y2, z2, color='darkgreen', linewidth=0.5)

# find intersections
for i in range(len(l2_coords)):
    if i == 0:  # skip first, there is no previous point
        continue

    distance, close_index = targ_kdtree.query(l2_coords[i], distance_upper_bound=.5)

    # strangely, points infinitely far away are somehow within the upper bound
    if np.isinf(distance):
        continue

    # plot ground truth that was activated
    _x, _y, _z = targ_kdtree.data[close_index]
    ax.scatter(_x, _y, _z, 'gx')
    _x2, _y2, _z2 = l2_coords[i]
    ax.scatter(_x2, _y2, _z2, 'bx')  # Plot the cross point


plt.show()