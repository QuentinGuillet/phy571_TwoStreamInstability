import sys
sys.path.append('../') 
import Core.Utility as util
#import Core.Savers as Savers
import numpy as np
import matplotlib.pyplot as plt

'''This program prints the position of a particule in the 2-symmetric and antisymmetric modes'''


'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 15000
dt = 0.001
n = 100
dx = 1./n
'''Initial conditions'''

N = 2
pos_init = np.array([0.25,0.75])
speed_init = np.array([.15,-.15])

plasma1 = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
for k in range(T):
   plasma1.move_one_turn()
   





import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10.0


'''Print the trajectories'''

plt.xlabel("Time (t*$\omega_p$)",fontsize="large")
plt.ylabel("Position (x/L)",rotation="vertical",position=(0,0.58))
print(plasma1.positionsstored[:,0])

ax, = plt.plot(np.linspace(0,T*dt,T),plasma1.positionsstored[:,0])
ax2, = plt.plot(np.linspace(0,T*dt,T),plasma1.positionsstored[:,1])
plt.legend((ax,ax2), ('Particle 1', 'Particle 2'),loc=1)


'''Print the trajectory in the phase space'''
'''
plt.xlabel("Position (x/L)",fontsize="large")
plt.ylabel("Speed (v/L$\omega_p$)",rotation="vertical",position=(0,0.58))
print(plasma1.positionsstored[:,0])

ax, = plt.plot(plasma1.positionsstored[:,0],plasma1.speedstored[:,0])
ax2, = plt.plot(plasma1.positionsstored[:,1],plasma1.speedstored[:,1])
plt.legend((ax,ax2), ('Particle 1', 'Particle 2'),loc=1)
'''


plt.show()
