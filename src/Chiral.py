import numpy as np
from alive_progress import alive_bar
import matrix_routines as Mat
import Exceptions
from Stats import Stats

class Chiral(object):
    
    def __init__(self, N, beta, N_measurment, N_thermal, N_sweeps, epsilon, N_tau, SU = 3, a = 1, order = 10, order_N = 10, renorm_freq = 1000, Hot_start = True, accel = False) -> None:
        
        self.SU = SU
        
        self.a = a
        self.accel = accel
        
        self.order_N = order_N
        
        self.renorm_freq = renorm_freq

        self.epsilon = epsilon
        
        self.N_tau = N_tau
        
        self.N_sweeps = N_sweeps
        self.mass = 0.1
        
        if SU == 3 or SU == 2 or SU==4:
        
            self.c = 2
        
        else:
        
            raise ValueError('Not considered value for SU')
        
        self.N = N
        ks = np.arange(0, self.N) # lattice sites in Fourier space along one direction
        A = np.zeros((self.N,self.N)) # inverse kernel computed at every site in Fourier space
        for k in range(self.N):
            for k_ in range(k,self.N):
                A[k,k_] = ( 4*np.sin(np.pi*ks[k]/self.N)**2 + 4*np.sin(np.pi*ks[k_]/self.N)**2 + self.mass**2)**(-1)   
                A[k_,k] = A[k,k_] # exploit symmetry of kernel under exchange of directions 
        self.A = A
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

    def Sample_FA_momentum(self):
        pi_F = np.zeros((self.N, self.N, self.SU**2-1), dtype=complex)

        PI_std = np.sqrt(self.N**2 / self.A) 
        STD = np.repeat(PI_std[:,:,None], repeats=self.SU**2-1, axis=2) # standard deviation is identical for components at same position
        PI = np.random.normal(loc=0, scale=STD) #  (N,N,self.SU**2-1) as returned array matches shape of STD

        # assign special modes for which FT exponential becomes +/-1. To get real pi in real space, the modes must be real themselves.
        N_2 = int(self.N/2)
        # two spacial indices
        pi_F[0,0] = PI[0,0]
        pi_F[0,N_2] = PI[0,N_2]
        pi_F[N_2,0] = PI[N_2,0]
        pi_F[N_2,N_2] = PI[N_2,N_2]

        # one special index
        pi_F[0,1:N_2] = 1/np.sqrt(2) * (PI[0,1:N_2] + 1j * PI[0,N_2+1:][::-1])
        pi_F[0,N_2+1:] = np.conj(pi_F[0,1:N_2][::-1]) # imposing hermitean symmetry

        pi_F[N_2,1:N_2] = 1/np.sqrt(2) * (PI[N_2,1:N_2] + 1j * PI[N_2,N_2+1:][::-1])
        pi_F[N_2,N_2+1:] = np.conj(pi_F[N_2,1:N_2][::-1])

        pi_F[1:N_2,0] = 1/np.sqrt(2) * (PI[1:N_2,0] + 1j * PI[N_2+1:,0][::-1])
        pi_F[N_2+1:,0] = np.conj(pi_F[1:N_2,0][::-1])

        pi_F[1:N_2,N_2] = 1/np.sqrt(2) * (PI[1:N_2,N_2] + 1j * PI[N_2+1:,N_2][::-1])
        pi_F[N_2+1:,N_2] = np.conj(pi_F[1:N_2,N_2][::-1])

        # no special index
        pi_F[1:N_2,1:N_2] = 1/np.sqrt(2) * (PI[1:N_2,1:N_2] + 1j * PI[N_2+1:,N_2+1:][::-1,::-1])
        pi_F[N_2+1:,N_2+1:] = np.conj(pi_F[1:N_2,1:N_2][::-1,::-1]) # imposing hermitean symmetry
   
        pi_F[1:N_2,N_2+1:] = 1/np.sqrt(2) * (PI[1:N_2,N_2+1:] + 1j * PI[N_2+1:,1:N_2][::-1,::-1])
        pi_F[N_2+1:,1:N_2] = np.conj(pi_F[1:N_2,N_2+1:][::-1,::-1])

        # pi is real by construction
        pi = np.real(np.fft.ifft2(pi_F, axes=(0,1)))

        return pi
    @staticmethod
    def Hamiltonian_FA(p,U,beta,c,A):
        N = len(U)
        p1 = Mat.dagger(np.fft.fft2( p,axes=(0,1)))
        p2  = np.einsum('ij,ijlk->ijlk', A, np.fft.fft2(p,axes=(0,1)))

        return 1/(2 *(c*N**2)) * np.einsum('ijkl,ijlk->', p1, p2).real+ Chiral.action(U, beta,c)
    @staticmethod
    def molecular_dynamics_FA(p_0, U_0, beta, epsilon, N_tau, SU, order, identity, order_N,A):
        p = p_0 + epsilon/2 * Chiral.dot_p(U_0,beta, SU, identity)
        def p_f(A,ps):
            mult = np.einsum('ij,ijlk->ijlk',A, np.fft.fft2(ps,axes=(0,1)))
            return np.fft.ifft2(mult,axes=(0,1))
        exponential_matrix = Mat.exponential(epsilon * p_f(A,p), order, SU, order_N)
        
        Unew =np.matmul(exponential_matrix, U_0)
        #np.einsum('ijkl,ijlm->ijkm', exponential_matrix,U_0)
       
        for i in range(N_tau):
            
            p += epsilon * Chiral.dot_p(Unew, beta, SU, identity)
            
            exponential_matrix = Mat.exponential(epsilon * p_f(A,p), order, SU, order_N)
            
            
            Unew = np.matmul(exponential_matrix, Unew)
            #np.einsum('ijkl,ijlm->ijkm', exponential_matrix, Unew)
        
        p += epsilon/2 * Chiral.dot_p(Unew, beta, SU, identity)
        
        return p, Unew

    def HMC_FA(self, flag = True):
        ## Sampling
        p_i = self.Sample_FA_momentum()
        p = np.einsum('abi,ikl->abkl', p_i, self.generators)

        ## Obtain Hamiltonian
        H = Chiral.Hamiltonian_FA(p, self.U, self.beta, self.c,self.A)

        ## Molecular Dynamics

        p_new, U_new = Chiral.molecular_dynamics_FA(p.copy(), self.U.copy(), self.beta, self.epsilon, self.N_tau, self.SU, self.order, self.identity, self.order_N,self.A)
        p_new = -p_new
        H_new = Chiral.Hamiltonian_FA(p_new, U_new, self.beta, self.c,self.A)
        
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
                if self.accel:
                    self.HMC_FA()
                else:
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

                    if self.accel:
                        self.HMC_FA()
                    else:
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
            print(np.average(np.exp(-self.DH)))

  
        print('------------------------------')
        
        return results, rate
    
    def Calibration_Runs(self, N_runs, N_thermal):
        """
        Runs the HMC for Calibration of the parameters of the algorithms
        """
        DH = np.zeros(N_runs)
    
        with alive_bar(N_runs + N_thermal) as bar:
            
            for i in range(N_thermal):
       
                if self.accel:
                        self.HMC_FA()
                else:
                    self.HMC()
           
                bar()
                
           
            self.accepted = 0
           
            self.tries = 0
            
            for i in range(N_runs):
                
                if self.accel:
                        self.HMC_FA()
                else:
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

