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
def zero(x,y):
    return 0 
    
class Plasma(object):
     
    def save_to_array(self):
        self.positionsstored[self.timestamp] = self.pos
        self.speedstored[self.timestamp] = self.speed 
        self.energystored[self.timestamp] = self.energy   
        #self.electric[self.timestamp] = self.E     
    
    
    def compute_energy(self):
        
        j = (self.pos/self.dx).astype(int)
        factor = (self.pos - self.X[[j[:,0]],[j[:,1]]][0])/4
        phi_part = (\
                    (1./2 - factor[:,0] - factor[:,1])*(self.phi + self.V)[[j[:,0]],[j[:,1]]]/self.dx + \
                    (1./4 - factor[:,0] + factor[:,1])*(self.phi + self.V)[[j[:,0]],[(j[:,1]+1)%self.n]]/self.dx + \
                    (1./4 + factor[:,0] - factor[:,1])*(self.phi + self.V)[[(j[:,0]+1)%self.n],[j[:,1]]]/self.dx + \
                    (1./4 + factor[:,0] - factor[:,1])*(self.phi + self.V)[[(j[:,0]+1)%self.n],[(j[:,1]+1)%self.n]]/self.dx)[0]
        
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
        
        self.X = np.mgrid[0:1:(n+1)*1j, 0:1:(n+1)*1j].transpose(1,2,0)
        self.V = np.vectorize(V)(self.X[:-1,:-1,0],self.X[:-1,:-1,1])
        self.positionsstored = np.zeros((T,N,2))
        self.speedstored = np.zeros((T,N,2))
        self.energystored = np.zeros(T)
        
        '''Creating the particules and information arrays'''
        self.pos = pos_init #pos_init should be a numpy array with dimension N (at time t=0)
        self.speed = speed_init #at time t = dt/2, see the leapfrog algorithm
        self.rho = np.zeros((n,n))
        self.rho_ = np.zeros((n,n),dtype=np.complex64)
        self.phi = np.zeros((n,n))
        self.phi_ = np.zeros((n,n),dtype=np.complex64)
        self.E = np.zeros((n,n))
        #self.compute_energy()
        #self.save_to_array()
        
    
    def position_to_density(self): #Transforms the particule positions to a density array along the grid
        self.rho = np.zeros((self.n,self.n))
        indices_in_grid = (self.pos/self.dx).astype(int)
        X_sorted = (self.X[[indices_in_grid[:,0]],[indices_in_grid[:,1]]])[0]
        rho_values = (self.pos-X_sorted)/self.dx
        #Affectation des coeeficients sur les 4 points entourant la particule
        np.add.at(self.rho,(indices_in_grid[:,0],indices_in_grid[:,1]),                         (1-rho_values[:,0])*(1-rho_values[:,1])/self.dx)
        np.add.at(self.rho,((indices_in_grid[:,0]+1)%self.n,indices_in_grid[:,1]),              (rho_values[:,0])*(1-rho_values[:,1])/self.dx)
        np.add.at(self.rho,(indices_in_grid[:,0],(indices_in_grid[:,1]+1)%self.n),              (1-rho_values[:,0])*(rho_values[:,1])/self.dx)
        np.add.at(self.rho,((indices_in_grid[:,0]+1)%self.n,(indices_in_grid[:,1]+1)%self.n),   (rho_values[:,0])*(rho_values[:,1])/self.dx)
        
    def compute_ElectricField(self):
        self.rho_ = np.fft.fft2(self.rho)
        #after this operation whe have rho_(k,l) = DFT(k,l)(self.rho), we have to compute the inverse of the laplacian accordingly
        # That is exactly design an array K[i,j] = (i**2 + j**2)*2*np.pi/self.n of shape (self.n,self.n) and divide our fourrier transform by this array
        # This K-array becomes obsolete then K = np.fft.fftfreq(self.n*self.n,self.dx**2).reshape(self.n,self.n)*2*np.pi
        K = np.mgrid[0:1-1./self.n:self.n*1j, 0:1-1./self.n:self.n*1j].transpose(1,2,0)*2*np.pi
        self.phi_[1:,1:] = self.rho_[1:,1:]/self.eps0/(np.linalg.norm(K,axis=1)[1:]**2)
        self.phi = np.fft.ifft2(self.phi_).real + self.V
        self.E = np.zeros((self.n,self.n,2))
        B = (np.roll(self.phi,1,0)-np.roll(self.phi,-1,0))/2/self.dx
        C = (np.roll(self.phi,1,1)-np.roll(self.phi,-1,1))/2/self.dx
        self.E[:,:,0] = B
        self.E[:,:,1] = C
   
    def field_to_particules(self): #Transforms a field on the grid to a force applied to particules, returns the array of forces on the particules
        Force = self.q*self.E
        j = (self.pos/self.dx).astype(int)
        factor = ((self.pos - self.X[[j[:,0]],[j[:,1]]])/4/self.dx)[0]
        multiply = np.zeros((self.N,2,2))
        multiply[:,0,0]  = factor[:,0]
        multiply[:,0,1]  = factor[:,0]
        multiply[:,1,0]  = factor[:,1]
        multiply[:,1,1]  = factor[:,1]
        
        # f ( x , y ) ≈ f ( 0 , 0 ) ( 1 − x ) ( 1 − y ) + f ( 1 , 0 ) x ( 1 − y ) + f ( 0 , 1 ) ( 1 − x ) y + f ( 1 , 1 ) x y 
        #This is bilinear interpolation
        self.F =    Force[[j[:,0]],[j[:,1]]][0]*(1-multiply[:,0])*(1-multiply[:,1]) + \
                    Force[[j[:,0]],[(j[:,1]+1)%self.n]][0]*(1-multiply[:,0])*(multiply[:,1]) + \
                    Force[[(j[:,0]+1)%self.n],[j[:,1]]][0]*(multiply[:,0])*(1-multiply[:,1]) + \
                    Force[[(j[:,0]+1)%self.n],[(j[:,1]+1)%self.n]][0]*(multiply[:,0])*(multiply[:,1])
                   
           
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
        np.save("../results/2Dposition.npy", self.positionsstored)
        np.save("../results/2Dspeed.npy", self.speedstored)
        np.save("../results/2Denergy.npy", self.energystored)
        #np.savetxt("../results/field.dat", self.electric)
