from ex2 import grav as gv

from matplotlib import pyplot as plt

gv.seed = 1369420
gv.ti = 0
gv.dt = 0.001

gv.main()

t = gv.t[:]
ec = gv.ec[:]
ep = gv.ep[:]
et = gv.et[:]

plt.figure()
plt.title('Variation of system energy over time')
plt.xlim(-100,1100)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.plot(t,ec,label='Kinetic Energy')
plt.plot(t,ep,label='Potential Energy')
plt.plot(t,et,label='Total Energy')
plt.legend(fontsize='small')
plt.savefig('ex2_energy_l100_p100.png',dpi=300)
plt.close()
