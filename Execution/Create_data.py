import sys
sys.path.append('../') 
import Core.Utility as util
import numpy as np

'''[Constants]'''
q = 1
m = 1
eps0 = 1
T = 1000
n = 10000
dt = 0.01
dx = 1/n




'''[Initial Conditions]'''
'''For two particules'''
'''
N = 2
pos_init = np.array([0.25,0.75])
speed_init = np.array([0.1,-0.1])
'''


'''For N particules equally repartitionned'''
N = 1000
pos_init_1 = np.linspace(1/(2*N),1-1/(2*N),N)
pos_init_2 = np.linspace(0,1-1/N,N)
speed_init_1 = np.zeros(N) + 0.1
pos_init = np.concatenate((pos_init_1,pos_init_2))
speed_init = np.concatenate((speed_init_1,-speed_init_1))
N*=2
print(pos_init)
plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
for i in range(T):
    plasma.move_one_turn()
plasma.save_to_text()
