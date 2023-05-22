from grav_2 import grav as gv

from matplotlib import pyplot as plt

from celluloid import Camera

import numpy as np

gv.seed = 1369420
gv.ti = 0
gv.dt = 0.001

n = int(gv.n)


gv.main()

t = gv.t[:]

fig = plt.figure()
ax = fig.gca(projection='3d')
camera = Camera(fig)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

X = np.zeros(int(gv.p))
Y= np.zeros(int(gv.p))
Z = np.zeros(int(gv.p))

k = 0
while k <= n:
    for i in range(len(X)): # Gravar a posição em j das p partículas.
        x = gv.x[0,:,i]
        y = gv.x[1,:,i]
        z = gv.x[2,:,i]
        X[i] = x[k-1]
        Y[i] = y[k-1]
        Z[i] = z[k-1]
    k += 10000
    ax.scatter(X,Y,Z,marker='o',c='black')
    camera.snap()
animation = camera.animate()
animation.save('ex2_animation_l10_p10.mp4')
