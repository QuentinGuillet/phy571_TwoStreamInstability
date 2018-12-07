import sys
sys.path.append('../') 
import Core.Utility as util
#import Core.Savers as Savers
import numpy as np
import matplotlib.pyplot as plt
from random import gauss



'''This program computes the energy error and plots it against n, the grid decomposition in the instability case, after the instability (N=10000)'''

'''The instability happens before t=0.4s, so we keep Tmax*dt = 0.4'''

'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 400
dt = 0.001

    
'''Gaussian distributions'''
initial_speed_mean = 10
initial_speed_sigma = 1

ndebut = 20
nfin = 200
spacing = 5

dE = np.zeros((nfin - ndebut)//spacing)



for i in range(ndebut+1,nfin,spacing):
    n = i
    print(i)
    dx = 1./n
    
    
    '''Initial Conditions'''
    N = 10000
    pos_init_1 = np.linspace(1./(2*N),1.-1./(2*N),N)
    pos_init_2 = np.linspace(0,1.-1./N,N)
    pos_init = np.concatenate((pos_init_1,pos_init_2))
    
    speed_init_1 = np.array([gauss(initial_speed_mean,initial_speed_sigma) for i in range(N)])
    speed_init_2 = np.array([gauss(-initial_speed_mean,initial_speed_sigma) for i in range(N)])
    speed_init = np.concatenate((speed_init_1,speed_init_2))
    N*=2
        
        

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
ax, = plt.plot(np.linspace(ndebut+1,nfin,(nfin - ndebut )//spacing - 1),dE[1:],"tomato")
plt.show()
