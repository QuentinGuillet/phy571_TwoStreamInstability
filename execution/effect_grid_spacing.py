import sys
sys.path.append('../') 
import Core.Utility as util
#import Core.Savers as Savers
import numpy as np
import matplotlib.pyplot as plt



'''This program computes the figure of final speed depending on n'''


'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 1500
dt = 0.001


nx = 30
vfin = np.zeros(nx)
for i in range(1,nx):
    n = i
    dx = 1./n
    N = 2
    pos_init = np.array([0.25,0.75])
    speed_init = np.array([0.0,-0.0])
    plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
    print(speed_init)
    for k in range(T):
       plasma.move_one_turn()
    vfin[i] = plasma.speedstored[-1][0]


import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10.0


plt.xlabel("Grid decomposition (n)",fontsize="large")
ax, = plt.plot(vfin)
plt.legend((ax,),("Max speed",))
plt.show()
