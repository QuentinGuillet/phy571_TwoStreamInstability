import numpy as np


q = -1   #Charge
deltax = 0.5 #Spatial step
deltat = 0.5 #Temporal step
n = 1000 #Number of points on the grid
N = 10000 #Number of particles
T = 100 #Number of temporal steps
X = np.linspace(0,n*deltax,n) #Grid array
m = 0.1 #Mass of a particule


def position_to_density(pos): #Transforms the particule positions to a density array along the grid
    rho = np.zeros(n)
    for i in range(N):
        j,alpha,beta = barycentre(pos[i]) 
        rho[j] += alpha/deltax
        rho[(j+1)%n] += beta/deltax
    return rho

def barycentre(pos_i): #Returns the contribution of a particule at pos_i to the density on the grid
    j = int(pos_i/deltax)
    alpha = (pos_i-X[j])/deltax
    return j,1-alpha,alpha
 


def field_to_particules(E,pos): #Transforms a field on the grid to a force aplied to particules, returns the array of forces on the particules
    F = np.zeros(N)
    for i in range(N):
        F[i] = get_force(E,pos[i])
    return F


def get_force(E,pos_i,X): #Returns the force applied on point qi on the grid
    F = q*E
    j = int(pos_i/deltax)
    alpha = (X[(j+1)%n] - pos_i)*F[j]/deltax +  (pos_i - X[j])*F[j+1]/deltax
    return alpha
 
 
def move_particules(F,v,x): #Moves every particule according to the force applied and the speeds
    v+= F/m*deltat
    x+= v*deltat