from ex1 import verlet

import matplotlib.pyplot as plt

import numpy as np

from mpl_toolkits.mplot3d import Axes3D

verlet.dt = 0.01
verlet.t0 = 0.
verlet.k = 1000
verlet.lb = 1.5
verlet.xi = [-100,100,-100]
verlet.vi = [1,1,1]


verlet.vl()

x = verlet.x[0,:]
y = verlet.x[1,:]
z = verlet.x[2,:]
t = verlet.t[:]
ec = verlet.ec[:]
ep = verlet.ep[:]
et = verlet.et[:]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.plot(x, y, z)
ax.set_title('Space Solution for ' r'$\lambda = 1.5$')
ax.view_init(elev=10., azim=-25)
plt.savefig('space_solution_lambda15.png',dpi=300)
plt.close()

plt.figure()
plt.title('Solution in the 'r'$x \times y$' ' plane for ' r'$\lambda = 1.5$')
plt.xlabel('X')
plt.ylabel('Y')
plt.plot(x,y)
plt.savefig('x_times_y_plane_lambda15.png',dpi=300)
plt.close()

plt.figure()
plt.title('Energy ' r'$\times$' ' Time for ' r'$\lambda = 1.5$' )
plt.xlabel('Time')
plt.ylabel('Energy')
plt.ylim(min(ep) - 5, max(ec) + 15)
plt.plot(t,ec,label='Kinetic Energy')
plt.plot(t,ep,label='Potential Energy')
plt.plot(t,et,label='Total Energy')
plt.legend(fontsize='small')
plt.savefig('energy_lambda15.png',dpi=300)
plt.close()

