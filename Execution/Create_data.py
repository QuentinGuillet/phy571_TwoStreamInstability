import sys
sys.path.append('../') 
import Core.Utility as util
import numpy as np

'''[Constants]'''
q = 1
m = 1
eps0 = 1
T = 1000
n = 1000
dt = 0.01
dx = 1e-3
N = 2



'''[Initial Conditions]'''
pos_init = np.array([0.25,0.75])
speed_init = np.array([0.1,-0.1])




plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
for i in range(T):
    plasma.move_one_turn()
plasma.save_to_text()
