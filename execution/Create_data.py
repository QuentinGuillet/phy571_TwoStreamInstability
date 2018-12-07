import sys
sys.path.append('../') 
import Core.Utility as util
#import Core.Savers as Savers
import numpy as np
from random import gauss

'''[Constants]'''

'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 2000
n = 100
dt = 0.001
dx = 1./n

'''Gaussian distributions'''
initial_speed_mean = 10
initial_speed_sigma = 1

'''Potential'''
omega = 1

def Vt(x):
    return 1/2*omega**2*(x-0.5)**2

'''[Initial Conditions]'''
'''For two particules'''
'''
N = 2
pos_init = np.array([0.25,0.75])
speed_init = np.array([0.1,-0.1])
'''

'''For One particules'''
'''
N = 1
pos_init = np.array([0.5])
speed_init = np.array([0.])
'''


'''For N particules equally repartitionned'''

N = 10000
pos_init_1 = np.linspace(1./(2*N),1.-1./(2*N),N)
pos_init_2 = np.linspace(0,1.-1./N,N)
pos_init = np.concatenate((pos_init_1,pos_init_2))


'''All the particules have same initial speed'''
'''
speed_init_1 = np.zeros(N) + 3
speed_init_2 = - speed_init_1
'''

'''Gaussian speed distribution'''

speed_init_1 = np.array([gauss(initial_speed_mean,initial_speed_sigma) for i in range(N)])
speed_init_2 = np.array([gauss(-initial_speed_mean,initial_speed_sigma) for i in range(N)])
speed_init = np.concatenate((speed_init_1,speed_init_2))
N*=2


'''
Position_saver = Savers.PositionS("../results/position.npy")
Speed_saver = Savers.SpeedS("../results/speed.npy")
Energy_saver = Savers.EnergyS("../results/energy.npy")
Potential_saver = Savers.PotentialS("../results/potential.npy")
'''


plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
for i in range(T):
    plasma.move_one_turn()
plasma.save_to_text()