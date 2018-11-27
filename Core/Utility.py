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
def zero(x):
    return 0 
    
class Plasma(object):
    
   
    def save_to_array(self):
        self.positionsstored[self.timestamp] = self.pos
        self.speedstored[self.timestamp] = self.speed 
        self.energystored[self.timestamp] = self.energy     
        #self.electric[self.timestamp] = self.E     
    
    def compute_energy(self):
        phi_part = np.zeros(self.N)
        for i in range(self.N):
            #We assume here that 0 < pos_i < self.dx*self.n so that 0<= j < n
            pos_i = self.pos[i]
            j = int(pos_i/self.dx)
            if j>=self.n:
                j=0
            phi_part[i] = (self.X[(j+1)] - pos_i)*(self.phi[j]+self.V[j])/self.dx +  (pos_i - self.X[j])*(self.phi[(j+1)%self.n]+self.V[(j+1)%self.n])/self.dx
        self.energy = np.sum(0.5*self.m*self.speed**2) + np.sum(0.5 * self.q*phi_part)
             
        
    def __init__(self,q,m,dx,dt,n,N,T,eps0,pos_init,speed_init,V=zero):
        '''Initialization of the intern variables'''
        self.q = q
        self.m = m
        self.dx = dx
        self.dt = dt
        self.n = n
        self.N = N
        self.T = T
        self.eps0 = eps0
        self.timestamp = 0
        '''Initialization of the array variables'''
        self.X = np.linspace(0,n*dx,n+1) #We take the last point as a point of X to be sure the force at the limits is not absurd
        self.V = np.vectorize(V)(self.X[:n])
        print(self.V)
        self.positionsstored = np.zeros((T,N))
        self.speedstored = np.zeros((T,N))
        self.energystored = np.zeros(T)
        self.electric = np.zeros((T,n))
        
        '''Creating the particules and information arrays'''
        self.pos = pos_init #pos_init should be a numpy array with dimension N (at time t=0)
        self.speed = speed_init #at time t = dt/2, see the leapfrog algorithm
        self.rho = np.zeros(n)
        self.rho_ = np.zeros(n,dtype=np.complex64)
        self.phi = np.zeros(n)
        self.phi_ = np.zeros(n,dtype=np.complex64)
        self.E = np.zeros(n)
        self.compute_energy()
        self.save_to_array()
        
        
    def barycentre(self,pos_i): #Returns the contribution of a particule at pos_i to the density on the grid
        #We assume here that 0 < pos_i < self.dx*self.n so that 0<= j < n
        j = int(pos_i/self.dx)
        if j>=self.n:
            j=0
        alpha = (pos_i-self.X[j])/self.dx
        return j,1-alpha,alpha
    
    def position_to_density(self): #Transforms the particule positions to a density array along the grid
        self.rho = np.zeros(self.n)
        j = (self.pos/self.dx).astype(int)
        Xscaled = self.X[j]
        j = np.concatenate((j,(j+1)%self.n))
        np.add.at(self.rho,j,np.concatenate(((1-(self.pos-Xscaled)/self.dx)/self.dx,(self.pos-Xscaled)/self.dx/self.dx)))
        '''
        for i in range(self.N):
            j,alpha,beta = self.barycentre(self.pos[i]) 
            self.rho[j] += alpha/self.dx
            self.rho[(j+1)%self.n] += beta/self.dx
        #self.rho -= 1
        '''
            
    def compute_ElectricField(self):
        self.rho_ = np.fft.fft(self.rho)
        K = np.fft.fftfreq(self.n,self.dx)*2*np.pi
        self.phi_[1:] = self.rho_[1:]/self.eps0/K[1:]**2
        self.phi = np.fft.ifft(self.phi_).real + self.V
        self.E = (np.roll(self.phi,1)-np.roll(self.phi,-1))/2/self.dx        # E(x) = - (phi(x+dx)-phi(x-dx)/2dx)
 
 
    def field_to_particules(self): #Transforms a field on the grid to a force applied to particules, returns the array of forces on the particules
        self.F = np.zeros(self.N)
        Force = self.q*self.E
        for i in range(self.N):
            #We assume here that 0 < pos_i < self.dx*self.n so that 0<= j < n
            pos_i = self.pos[i]
            j = int(pos_i/self.dx)
            if j>=self.n:
                j=0
            self.F[i] = (self.X[(j+1)] - pos_i)*Force[j]/self.dx +  (pos_i - self.X[j])*Force[(j+1)%self.n]/self.dx
    
    
    def move_particules(self): #Moves every particule according to the new speed
        self.pos += self.speed*self.dt
        self.pos %= (self.dx*self.n)
        
        
    def compute_speed_and_energy(self): #This function is a bit weird because we need v(t) and x(t) simultaenously for a a good computation then we compute v(t+dt/2)
                                        #This is due to the leapfrog method
        self.speed += self.F/self.m*self.dt/2
        self.compute_energy()
        self.speed += self.F/self.m*self.dt/2
                        
    def move_one_turn(self): #Executes a loop of the principle
        if(self.timestamp < self.T - 1):
            self.move_particules()      # we move from to t to t + dt using v(t+dt/2)
            self.position_to_density()
            self.compute_ElectricField()
            self.field_to_particules()  
            self.compute_speed_and_energy()    # we compute v(t+dt/2+ dt) using x(t+dt)
            self.timestamp += 1
            self.save_to_array()
    
        
    def save_to_text(self):
        np.save("../results/position.npy", self.positionsstored)
        np.save("../results/speed.npy", self.speedstored)
        np.save("../results/energy.npy", self.energystored)
        #np.savetxt("../results/field.dat", self.electric)
