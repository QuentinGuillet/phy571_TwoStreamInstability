import sys
sys.path.append('../') 
import Core.Utility as util
#import Core.Savers as Savers
import numpy as np
import matplotlib.pyplot as plt



'''This program computes the energy error and plots it against n, the grid decomposition, in the 2 particle case'''


'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 1500
dt = 0.001

ndebut = 20
nfin = 200
spacing = 5
dE = np.zeros((nfin - ndebut)//spacing)
for i in range(ndebut+1,nfin,spacing):
    n = i
    print(i)
    dx = 1./n
    N = 2
    pos_init = np.array([0.25,0.75])
    speed_init = np.array([0.1,-0.1])
    plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
    for k in range(T):
       plasma.move_one_turn()
    energy = plasma.energystored[1:]
    dE[(i - ndebut - 1)//spacing] = float(max(energy) - min(energy))/float(sum(energy))*float(len(energy))

np.save("../results/energy_error.npy",dE)


import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10.0


plt.xlabel("Grid decomposition (n)",fontsize="large")
plt.ylabel("$\Delta E / E $",rotation="vertical",position=(0,0.58))
plt.plot(np.linspace(ndebut+1,nfin,(nfin - ndebut )//spacing - 1),dE[1:])
plt.show()
