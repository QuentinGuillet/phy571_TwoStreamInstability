'''This function returns the index of the grid point just below xi and the contributions to p(Xj) and p(Xj+1) for the i particule'''
import numpy as np
q = -1   #Charge
deltax = 0.5 #Spatial step
n = 1000 #Number of points on the grid
N = 10000 #Number of particles


def position_to_density(x):
    rho = np.zeros(n)
    for i in range(N):
        j,alpha,beta = barycentre(x[i]) 
        rho[j] += alpha
        rho[(j+1)%n] += beta
    return rho



def barycentre(pos_i):
 j = int(pos_i/deltax)
 alpha = (pos_i-j*deltax)/deltax
 return j,1-alpha,alpha
 
print(barycentre(0.4))

def get_force(E,qi,X):
    #   F = force_from_field(E)
 F = q*E
 j = int(qi/deltax)
 alpha = (X[j+1] - qi)*F[j]/deltax +  (qi - X[j])*F[j+1]/deltax
 return alpha