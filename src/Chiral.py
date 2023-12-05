import numpy as np
from alive_progress import alive_bar
import matrix_routines as Mat
import Exceptions
from Stats import Stats

class Chiral(object):
    
    def __init__(self, N, beta, N_measurment, N_thermal, N_sweeps, epsilon, N_tau, SU = 3, a = 1, order = 10, order_N = 10, renorm_freq = 1000, Hot_start = True) -> None:
        
        self.SU = SU
        
        self.a = a
        
        self.order_N = order_N
        
        self.renorm_freq = renorm_freq

        self.epsilon = epsilon
        
        self.N_tau = N_tau
        
        self.N_sweeps = N_sweeps
        
        if SU == 3 or SU == 2 or SU==4:
        
            self.c = 2
        
        else:
        
            raise ValueError('Not considered value for SU')
        
        self.N = N
        
        self.N_measurement = N_measurment
        
        self.N_thermal = N_thermal
        
        self.U = np.zeros((N, N, SU, SU), dtype = complex) 
        
        self.beta = beta
        
        self.generators = Mat.create_generators(self.SU)
        
        self.DH = np.zeros(self.N_measurement)
        
        self.order = order
        
        for i in range(N):
            
            for j in range(N):
            
                for k in range(SU):
            
                    self.U[i, j, k, k]  = 1
        
        self.identity = self.U.copy()
        
        self.accepted = 0
        
        self.tries = 0
        ## set-up Hot Start
        self.Hot_start = Hot_start
        if Hot_start:
            for i in range(10):

                self.HMC(False)

        self.accepted = 0
        
        self.tries = 0

        if N_thermal != 0:

            self.thermalize()
        
        self.accepted = 0
        
        self.tries = 0
    
    @staticmethod
    def action(U, beta, c):

        
        top_neighbours = np.roll(U, -1, axis = 0)
        
        right_neighbours = np.roll(U, 1, axis = 1)
        
        tops = Mat.multiply_matrices(Mat.dagger(U), top_neighbours)

        right = Mat.multiply_matrices(Mat.dagger(U), right_neighbours)

        return -beta/c * 2 * np.einsum('ijkk->',right + tops).real
    
   
    @staticmethod
    def Hamiltonian(p, U, beta, c):
        
        return 1/(2 * c) * np.einsum('ijkl,ijlk->', p, p).real+ Chiral.action(U, beta, c)

    @staticmethod
    def exponential_matrix(matrix, order):
        
        return np.sum([np.linalg.matrix_power(matrix, n)/np.math.factorial((n)) for n in range(order)], axis = 0, dtype = matrix.dtype)
    
    @staticmethod
    def dot_p(U, beta, SU, identity):
    
        top_neighbours = np.roll(U, -1, axis = 0)
    
        bottom_neighbours = np.roll(U, 1, axis = 0)
    
        left_nighbours = np.roll(U, -1, axis = 1)
    
        right_neighbours = np.roll(U, 1, axis = 1)
    
        vertical = Mat.multiply_matrices(top_neighbours + bottom_neighbours, Mat.dagger(U))
    
        horizontal = Mat.multiply_matrices(left_nighbours + right_neighbours, Mat.dagger(U))
    
        result = -1j * beta * (vertical -Mat.dagger(vertical) + horizontal - Mat.dagger(horizontal))
    
        if SU == 2:
    
            return result
    
        else:
    
            result-= np.einsum('ijkl,ijmm->ijkl', identity/SU, result)
    
        return result
    
    @staticmethod
    def molecular_dynamics(p_0, U_0, beta, epsilon, N_tau, SU, order, identity, order_N):
        
        p = p_0 + epsilon/2 * Chiral.dot_p(U_0,beta, SU, identity)
        
        exponential_matrix = Mat.exponential(epsilon * p, order, SU, order_N)
        
        Unew =np.matmul(exponential_matrix, U_0)
        #np.einsum('ijkl,ijlm->ijkm', exponential_matrix,U_0)
       
        for i in range(N_tau):
            
            p += epsilon * Chiral.dot_p(Unew, beta, SU, identity)
            
            exponential_matrix = Mat.exponential(epsilon * p, order, SU, order_N)
            
            
            Unew = np.matmul(exponential_matrix, Unew)
            #np.einsum('ijkl,ijlm->ijkm', exponential_matrix, Unew)
        
        p += epsilon/2 * Chiral.dot_p(Unew, beta, SU, identity)
        
        return p, Unew
    
    def HMC(self,flag=True):
        
        p_i = np.random.standard_normal(size = (self.N, self.N, int(self.SU**2 - 1)))
    
        p = np.einsum('abi,ikl->abkl', p_i, self.generators)

        H = Chiral.Hamiltonian(p, self.U, self.beta, self.c)
        
        p_new, U_new = Chiral.molecular_dynamics(p.copy(), self.U.copy(), self.beta, self.epsilon, self.N_tau, self.SU, self.order, self.identity, self.order_N)
        
        H_new = Chiral.Hamiltonian(p_new, U_new, self.beta, self.c)
        
        Delta_H = H_new - H
        
        self.delta_H = Delta_H
        if flag:
            if Delta_H <0 or np.exp(-Delta_H) > np.random.random():
        
                self.accepted += 1
            
                self.U = U_new.copy()
        else:
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
        if np.all(self.U == self.identity):
            raise Exceptions.ThermalizationException("The Field is still equal to the Identity when initialisating the Measurements, thermalisation has not occurred, possible fixes include running the program again, calibrating or increasing the number of thermalisation runs")
        results = [0 for i in range(self.N_measurement)]
        
        print('Measurements with beta = ' + str(self.beta) + " N = " +str(self.N))
        
        with alive_bar(self.N_measurement) as bar:


            for i in range(self.N_measurement):

                for j in range(self.N_sweeps):

                    self.HMC()
                
                self.DH[i]  = self.delta_H 
                
                results[i] = observable(self.U)
               
                bar()
                if i == 999 and self.accepted/self.tries <= 0.75:
                    raise Exceptions.CalibrationException('The Acceptance rate of the run is too low to be acceptable, consider recalibrating or running again')
                if (i%self.renorm_freq==0 and (self.SU == 3 or self.SU==4)):
                    self.U = Mat.reunitarisation(self.U.copy(), self.SU)
                    #print(np.average(Mat.determinant(self.U)))
            rate = self.accepted/self.tries
            
            print(rate)
  
        print('------------------------------')
        
        return results, rate
    
    def Calibration_Runs(self, N_runs, N_thermal):
        """
        Runs the HMC for Calibration of the parameters of the algorithms
        """
        DH = np.zeros(N_runs)
    
        with alive_bar(N_runs + N_thermal) as bar:
            
            for i in range(N_thermal):
       
                self.HMC()
           
                bar()
                
           
            self.accepted = 0
           
            self.tries = 0
            
            for i in range(N_runs):
                
                self.HMC()
                
                DH[i] = self.delta_H
                
                bar()
                
            
            print((self.accepted/self.tries)*100)
                    
    @staticmethod
    def Measure_Sus(U,SU):
        N = len(U)
        return 1/(SU*N**2) *(np.einsum('ijkl,mnlk->',U,Mat.dagger(U))).real
    @staticmethod
    def Measure_G(U,SU):
        Udag = Mat.dagger(U)
        N = len(U)
        result = np.zeros((N,N))
        result2 = np.zeros((N,N))
        for i in range(SU):
            for j in range(SU):
                A = np.fft.fft2(U[:,:,i,j])
                B = np.fft.fft2(np.conjugate((Udag[:,:,j,i])))
                result += np.real(np.fft.ifft2(A*np.conjugate(B)))
       
        #return 1/(SU) *np.einsum('ijkl,lk->ij',U,Mat.dagger(U)[0,0]).real
        return 1/(SU*N**2) *result
    @staticmethod
    def Measure_G_0_mom(U,SU):
        Udag = Mat.dagger(U)
        N = len(U)
        result = np.zeros((N,N))
        for i in range(SU):
            for j in range(SU):
                A = np.fft.fft2(U[:,:,i,j])
                B = np.fft.fft2(np.conjugate((Udag[:,:,j,i])))
                result += np.real(np.fft.ifft2(A*np.conjugate(B)))
        result = result/(SU*N**2)
        return np.einsum('ij->j',result)
    @staticmethod
    def Measure_G_diags(U,SU):
        Udag = Mat.dagger(U)
        N = len(U)
        result = np.zeros((N,N))
        for i in range(SU):
            for j in range(SU):
                A = np.fft.fft2(U[:,:,i,j])
                B = np.fft.fft2(np.conjugate((Udag[:,:,j,i])))
                result += np.real(np.fft.ifft2(A*np.conjugate(B)))
        result = result/(SU*N**2)
        Greens_diag = np.zeros((N))
        for i in range(N):
            Greens_diag[i] =sum([result[(i-j)%N,j] for j in range(N)])
        return Greens_diag

    @staticmethod
    def Measure_G2(U,SU):
        return 1/SU *np.einsum('ijkl,mnlk->ijmn',U,Mat.dagger(U)).real
    
    @staticmethod
    def Measure_ww_corr(U,SU):
        N = len(U)
        dagU = Mat.dagger(U)
        us = np.sum(U,axis = 1)
        us2 = np.sum(dagU,axis = 1)
        ww_cor = np.zeros(N)
        for l in range(SU):
            for k in range(SU):
                cf = Stats.correlator(us[:,l,k], np.conjugate(us2[:,k,l]))
                ww_cor += cf.real
        ww_cor = ww_cor/N**2
        """result = np.zeros(N)
        for t in range(N):
            for t0 in range(N):
                result[t] += np.einsum('jkl,mlk->',U[t0],dagU[(t+t0)%N]).real
        result = result/N**2"""
        return ww_cor

