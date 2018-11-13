import numpy as np


'''q = -1   #Charge
dx = 0.5 #Spatial step
dt = 0.5 #Temporal step
n = 1000 #Number of points on the grid
N = 10000 #Number of particles
T = 100 #Number of temporal steps
X = np.linspace(0,n*dx,n) #Grid array
m = 0.1 #Mass of a particule'''


#Classe de gestion de l'algorithme
    
class Plasma:
    
    def __init__(self,q,m,dx,dk,dt,n,N,T,eps0,pos_init,speed_init):
        '''Initialization of the intern variables'''
        self.q = q
        self.m = m
        self.dx = dx
        self.dt = dt
        self.dk = dk
        self.n = n
        self.N = N
        self.T = T
        self.eps0 = eps0
        
        '''Initialization of the array variables'''
        self.X = np.linspace(0,n*dx,n)
        self.K = np.linspace(dk,(n+1)*dk,n)
        
        '''Creating the particules and information arrays'''
        self.pos = pos_init #pos_init should be a numpy array with dimension N
        self.speed = speed_init
        self.rho = np.zeros(n)
        self.rho_ = np.zeros(n)
        self.phi = np.zeros(n)
        self.phi_ = np.zeros(n)
          
        def barycentre(self,pos_i): #Returns the contribution of a particule at pos_i to the density on the grid
            j = int(pos_i/dx)
            alpha = (pos_i-self.X[j])/dx
            return j,1-alpha,alpha
        
        def position_to_density(self): #Transforms the particule positions to a density array along the grid
            self.rho = np.zeros(n)
            for i in range(N):
                j,alpha,beta = self.barycentre(self.pos[i]) 
                self.rho[j] += alpha/dx
                self.rho[(j+1)%n] += beta/dx
            
        
        def compute_ElectricField(self): #Computes it from the density
            self.rho_ = np.fft.fft(self.rho)
            self.phi_ = self.rho_/eps0/self.K**2
            self.phi = np.fft.ifft(self.phi_)
            self.E = (np.roll(self.phi,1)-np.roll(self.phi,-1))/2/self.dx        # E(x) = - (phi(x+dx)-phi(x-dx)/2dx)
          
         
        def get_force(self,pos_i): #Returns the force applied on point qi on the grid
            F = q*self.E
            j = int(pos_i/dx)
            alpha = (self.X[(j+1)%n] - pos_i)*F[j]/dx +  (pos_i - self.X[j])*F[j+1]/dx
            return alpha
        
        
        def field_to_particules(self): #Transforms a field on the grid to a force aplied to particules, returns the array of forces on the particules
            self.F = np.zeros(N)
            for i in range(N):
                self.F[i] = get_force(self.E,self.pos[i])
        
        
        def move_particules(self): #Moves every particule according to the force applied and the speeds
            self.speed += self.F/self.m*self.dt
            self.pos += self.speed*self.dt
            
        def move_one_turn(self):
            self.position_to_density()
            self.compute_ElectricField()
            self.field_to_particules()
            self.move_particules()
            
        def save_to_text(self):
            np.savetxt("results/position.dat", self.pos)
            np.savetxt("results/speed.dat", self.speed)