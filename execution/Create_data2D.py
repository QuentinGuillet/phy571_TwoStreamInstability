import sys
sys.path.append('../') 
import Core.Utility2D as util
import numpy as np
from random import gauss
import matplotlib.pyplot as plt

'''[Constants]'''

'''Physical constants'''
q = -1
m = 1
eps0 = 1

'''Simulation Constants'''
T = 200
n = 40
dt = 0.01
dx = 1./n

'''Gaussian distributions'''
Temp = 300
initial_speed_mean = 0.1
initial_speed_sigma = 0.001

'''Potential'''
omega = 1

def Vt(x,y):
    return 1/2*omega**2*(x-0.5)*(y-0.5)**2

'''[Initial Conditions]'''
'''For two particules'''
'''
N = 4
pos_init = np.array([[0.75,0.5],[0.25,0.5],[0.5,0.75],[0.5,0.25]])
speed_init = np.array([[0.0,0.],[-0.0,0.],[0.0,0.],[-0.0,0.]])
'''

'''For One particules'''

N = 2
pos_init = np.array([[0.5,0.40],[0.75,0.40]])
speed_init = np.array([[0.,0.],[0.,0.]])


'''For N particules equally repartitionned'''
'''
N = 10

pos_init_1 = np.linspace(1./(2*N),1.-1./(2*N),N)
pos_init = np.zeros((2*N,2))
pos_init[:N,0] = 0.33
pos_init[:N,1] = pos_init_1
pos_init[N:,0] = 0.667
pos_init[N:,1] = pos_init_1
'''
'''All the particules have same initial speed'''
'''
speed_init_1 = np.zeros((N,2)) + np.array([0.1,0])
speed_init_2 = - speed_init_1
'''
'''Gaussian speed distribution'''
'''
speed_init_1 = np.array([[gauss(initial_speed_mean,initial_speed_sigma),gauss(initial_speed_mean,initial_speed_sigma)] for i in range(N)])
speed_init_2 = np.array([[gauss(-initial_speed_mean,initial_speed_sigma),gauss(-initial_speed_mean,initial_speed_sigma)] for i in range(N)])
'''


'''
speed_init = np.concatenate((speed_init_1,speed_init_2))

N=2*N

'''


plasma = util.Plasma(q,m,dx,dt,n,N,T,eps0,pos_init,speed_init)
for i in range(T):
    plasma.move_one_turn()
plasma.save_to_text()
'''
z = plasma.potentialcurve.transpose()
fig = plt.figure(figsize=plt.figaspect(2.))
ax = fig.add_subplot(2, 1, 2, projection='3d')
x = np.linspace(0,plasma.T,plasma.T)
y = plasma.X[0,:-1,1]
X,Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y, z, rstride=1, cstride=1, linewidth=0, antialiased=False)
ax.set_zlim3d(-1,1)
plt.show()
'''