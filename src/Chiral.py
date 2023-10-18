import numpy as np
from alive_progress import alive_bar


class Chiral(object):
    def __init__(self,N, beta, N_measurment, N_thermal, N_sweeps, SU=3, a = 1) -> None:
        self.SU = 3
        self.N_sweeps = N_sweeps
        if SU == 3 or SU == 2:
            self.c = 2
        else:
            raise ValueError('Not considered value for SU')
        self.N = N
        self.N_measurement = N_measurment
        self.N_thermal = N_thermal
        self.U = np.zeros((N,N,2,2),dtype=float) 
        self.beta = beta
        self.generators = self.create_generators()
        self.DH = np.zeros(self.N_measurement)
        #print(generators)
        """for i in range(N):

            for j in range(N):
                    self.U[i,j] = [1,0,1,0]"""
        
        for i in range(N):
            for j in range(N):
                for k in range(2):
                    self.U[i,j,k,k]  = 1
        self.thermalize()

    def create_generators(self):
        generators = np.zeros((8,3,3),dtype=complex)
        generators[0][0,1] = 1
        generators[0][1,0] = 1
        generators[1][0,1] = -1j
        generators[1][1,0] = 1j
        generators[2][0,0]  = 1
        generators[2][1,1] = -1
        generators[3][0,2] = 1
        generators[3][2,0] = 1
        generators[4][0,2] = -1j
        generators[4][2,0] = 1j
        generators[5][1,2] = 1
        generators[5][2,1] = 1
        generators[6][1,2] = -1j
        generators[6][2,1] = 1j
        generators[7][0,0] = 1/np.sqrt(3)
        generators[7][1,1] = 1/np.sqrt(3)
        generators[7][2,2] = -2/np.sqrt(3)
        return generators

    

    @staticmethod
    def action(U, beta):
        G = 0
        top_neighbours = np.roll(U,-1,axis = 0)
        right_neighbours = np.roll(U,1,axis = 1)
        for t in range(len(U)):
            for x in range(len(U)):
                G+= np.dot(np.matrix(U[t,x]).getH(),top_neighbours[t,x])  +np.matrix(np.dot(np.matrix(U[t,x]).getH(),top_neighbours[t,x])).getH()
                G+= np.dot(np.matrix(U[t,x]).getH(),right_neighbours[t,x])  + np.matrix(np.dot(np.matrix(U[t,x]).getH(),right_neighbours[t,x])).getH()
        return -beta/2* np.trace(G)
    @staticmethod
    def Hamiltonian(p,U,beta, c):
        return 1/(2*c)*np.trace(np.matmul(p[:,:],p[:,:]),axis1=2,axis2=3) + Chiral.action(U,beta)
    @staticmethod
    def exponential_matrix(matrix, order):
        return np.sum([np.linalg.matrix_power(matrix,n)/np.math.factorial((n)) for n in range(order)],axis = 0,dtype = matrix.dtype)
    
    def dot_p(U, beta):
        G = np.zeros((len(U),len(U),3,3),dtype=complex)
        top_neighbours = np.roll(U,-1,axis = 0)
        bottom_neighbours = np.roll(U,1,axis = 0)
        left_nighbours = np.roll(U,-1,axis = 1)
        right_neighbours = np.roll(U,1,axis = 1)
        for t in range(len(U)):
            for x in range(len(U)):
                upper_a = np.dot(top_neighbours[t,x] + bottom_neighbours[t,x],np.matrix(U[t,x]).getH())
                sider_a =  np.dot(left_nighbours[t,x] + right_neighbours[t,x],np.matrix(U[t,x]).getH())
                G[t,x]= upper_a + np.matrix(upper_a).getH() + sider_a + np.matrix(sider_a).getH() 
        
        return -1j-beta*G
    
    def molecular_dynamics(p_0, U_0, beta, epsilon, N_tau):
        p = 0
        U = 0
        for i in range(N_tau):
            p = 0
            U = 0
        p = 0
        return p, U
    
    def HMC(self):
        p_i = np.random.normal(loc=0,scale=1,size=(self.N,self.N, 8))
        p_i = np.random.normal(loc=0,scale=1,size=(N,N, 8))
        p = np.zeros((self.N,self.N,self.SU,self.SU),dtype=complex)
        for i in range(self.N):
            for j in range(self.N):
                p[i,j] = sum(p_i[i,j,k] * self.generators[k] for k in range(int(self.SU**2-1)))
        
        H = Chiral.Hamiltonian(p,self.U,self.beta,self.c)
        p_new, U_new = Chiral.molecular_dynamics(p.copy(),self.U.copy(),self.beta,self.epsilon,self.N_tau)
        H_new = Chiral.Hamiltonian(p_new,U_new, self.beta, self.c)
        Delta_H = H_new - H
        self.delta_H = Delta_H
        if Delta_H <0 or np.exp(-Delta_H) > np.random.uniform(0,1):
       
            self.accepted += 1
        
            self.U = U_new.copy()
        self.tries += 1

    def thermalize(self):
        """
        Runs the HMC alogrithm for N_thermal times
        """
        with alive_bar(self.N_thermal) as bar:
        
            for i in range(self.N_thermal):
                self.HMC()
                bar()
    
    def generate_measurements(self, observable):
        """
        Runs the HMC N_measurment times and records the observable value in each measurement.
        """

        results = [0 for i in range(self.N_measurement)]
        
        print('Measurements with a = ' + str(self.a) + " N = " +str(self.N))
        
        with alive_bar(self.N_measurement) as bar:


            for i in range(self.N_measurement):

                for j in range(self.N_sweeps):

                    self.HMC()
                
                self.DH[i]  = self.delta_H 
                
                results[i] = observable(self.U)
               
                bar()
        
        print('------------------------------')
        
        return results




"""a = Chiral(2,2)
#print(Chiral.action(a.U,a.beta))
print(a.action(a.U,1))
#print(Chiral.action(a.U,a.beta))
"""
"""N = 100
generators = np.zeros((8,3,3),dtype=complex)
generators[0][0,1] = 1
generators[0][1,0] = 1
generators[1][0,1] = -1j
generators[1][1,0] = 1j
generators[2][0,0]  = 1
generators[2][1,1] = -1
generators[3][0,2] = 1
generators[3][2,0] = 1
generators[4][0,2] = -1j
generators[4][2,0] = 1j
generators[5][1,2] = 1
generators[5][2,1] = 1
generators[6][1,2] = -1j
generators[6][2,1] = 1j
generators[7][0,0] = 1/np.sqrt(3)
generators[7][1,1] = 1/np.sqrt(3)
generators[7][2,2] = -2/np.sqrt(3)
p_i = np.random.normal(loc=0,scale=1,size=(N,N, 8))
p = np.zeros((N,N,3,3))
for i in range(N):
    for j in range(N):
        p[i,j] = sum(p_i[i,j,k] * generators[k] for k in range(8))
"""