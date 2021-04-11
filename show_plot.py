from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
x = []
y = []
z = []

for line in open('test.dat', 'r'):
    values = [float(s) for s in line.split()]
    x.append(values[0])
    y.append(values[1])
    z.append(values[2])
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
plt.show()

