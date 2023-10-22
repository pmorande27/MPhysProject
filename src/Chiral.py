import numpy as np
from alive_progress import alive_bar
import Matrix_Routines as Mat


class Chiral(object):
    def __init__(self,N, beta, N_measurment, N_thermal, N_sweeps, epsilon, N_tau, SU=3, a = 1, order = 10) -> None:
        self.SU = SU
        self.a = a
        
        self.epsilon = epsilon
        self.N_tau = N_tau
        self.N_sweeps = N_sweeps

        self.randoms = np.load('Test 2.npy')
        self.ps = np.load('Test.npy')
        if SU == 3 or SU == 2:
            self.c = 2
        else:
            raise ValueError('Not considered value for SU')
        self.N = N
        self.N_measurement = N_measurment
        self.N_thermal = N_thermal
        self.U = np.zeros((N,N,SU,SU),dtype=complex) 
        self.beta = beta
        self.generators = Mat.create_generators(self.SU)
        self.DH = np.zeros(self.N_measurement)
        self.order = order
        #print(generators)
        """for i in range(N):

            for j in range(N):
                    self.U[i,j] = [1,0,1,0]"""
        
        for i in range(N):
            for j in range(N):
                for k in range(SU):
                    self.U[i,j,k,k]  = 1
        self.identity = self.U.copy()
        self.accepted = 0
        self.tries = 0
        if N_thermal != 0:

            self.thermalize()
        self.accepted = 0
        self.tries = 0
    



    

    @staticmethod
    def action(U, beta):
        top_neighbours = np.roll(U,-1,axis = 0)
        right_neighbours = np.roll(U,1,axis = 1)
        tops = Mat.multiply_matrices(Mat.dagger(U),top_neighbours)
        #tops += Mat.dagger(tops)
        right = Mat.multiply_matrices(Mat.dagger(U),right_neighbours)
        #right += Mat.dagger(right)
        return -beta *np.einsum('ijkk->',right+tops).real
    
   
    @staticmethod
    def Hamiltonian(p,U,beta, c):
        return 1/(2*c)*np.einsum('ijkl,ijlk->',p,p).real+ Chiral.action(U,beta)

    @staticmethod
    def exponential_matrix(matrix, order):
        return np.sum([np.linalg.matrix_power(matrix,n)/np.math.factorial((n)) for n in range(order)],axis = 0,dtype = matrix.dtype)
    @staticmethod
    def dot_p(U, beta):
        top_neighbours = np.roll(U,-1,axis = 0)
        bottom_neighbours = np.roll(U,1,axis = 0)
        left_nighbours = np.roll(U,-1,axis = 1)
        right_neighbours = np.roll(U,1,axis = 1)
        vertical = Mat.multiply_matrices(top_neighbours+bottom_neighbours,Mat.dagger(U))
        horizontal =  Mat.multiply_matrices(left_nighbours+right_neighbours,Mat.dagger(U))
        return -1j*beta*(vertical -Mat.dagger(vertical) + horizontal - Mat.dagger(horizontal))

    @staticmethod
    def molecular_dynamics(p_0, U_0, beta, epsilon, N_tau,SU, order, identity):
        p = p_0 +epsilon/2*Chiral.dot_p(U_0,beta)
        exponential_matrix = Mat.exponential(epsilon*p,order,SU,identity,10)
        Unew = np.einsum('ijkl,ijlm->ijkm',exponential_matrix,U_0)
       
        for i in range(N_tau):
            p += epsilon*Chiral.dot_p(Unew,beta)
            exponential_matrix = Mat.exponential(epsilon*p, order,SU,identity,10)
            Unew = np.einsum('ijkl,ijlm->ijkm',exponential_matrix,Unew)
        p += +epsilon/2*Chiral.dot_p(Unew,beta)
        return p, Unew
    
    def HMC(self):
        p_i = np.random.standard_normal(size=(self.N,self.N, int(self.SU**2-1)))
        #p_i = self.ps[i]
        #p = np.zeros((self.N,self.N,self.SU,self.SU),dtype=complex)
    
        p = np.einsum('abi,ikl->abkl',p_i,self.generators)

        H = Chiral.Hamiltonian(p,self.U,self.beta,self.c)
        p_new, U_new = Chiral.molecular_dynamics(p.copy(),self.U.copy(),self.beta,self.epsilon,self.N_tau, self.SU, self.order,self.identity)
        H_new = Chiral.Hamiltonian(p_new,U_new, self.beta, self.c)
        Delta_H = H_new - H
        self.delta_H = Delta_H
        if Delta_H <0 or np.exp(-Delta_H) > np.random.random():
       
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
        
        print('Measurements with beta = ' + str(self.beta) + " N = " +str(self.N))
        
        with alive_bar(self.N_measurement) as bar:


            for i in range(self.N_measurement):

                for j in range(self.N_sweeps):

                    self.HMC()
                
                self.DH[i]  = self.delta_H 
                
                results[i] = observable(self.U)
               
                bar()
            #print(np.average(np.exp(-self.DH)))
            print(self.accepted/self.tries)
            

        
        print('------------------------------')
        
        return results
    def Calibration_Runs(self, N_runs, N_thermal):
        """
        Runs the HMC for Calibration of the parameters of the algorithms
        """
        DH = np.zeros(N_runs)
        with alive_bar(N_runs+N_thermal) as bar:
            
            for i in range(N_thermal):
               
           
                self.HMC()
           
                bar()
           
            self.accepted = 0
           
            self.tries = 0
            
            for i in range(N_runs):
                self.HMC()
                DH[i] = self.delta_H
                bar()
            #print(np.average(np.exp(-DH)))
            print((self.accepted/self.tries)*100)
                    

