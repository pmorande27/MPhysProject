import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats
from scipy.optimize import curve_fit
from alive_progress import alive_bar


def analytical(a,N,mu,m):

    R = 1+ a**2*mu**2/(2*m)-a*mu/np.sqrt(m)*(1+a**2*mu**2/(4*m))**0.5

    return 1/(2*mu*(m+a**2*mu**2/4)**0.5)*(1+R**N)/(1-R**N)

class Lattice(object):
    """
    Class used to represent the 1 Dimensional Lattice of the Quantum Harmonic Oscillator
    """

    def __init__(self, N, a, N_thermal, N_measurement, N_sweeps, Delta, N_tau, d_tau, mass, w, acceleration = False) -> None:
        """
        Constructor of the Class
        """
        self.delta_H = 0
        self.m = mass 

        self.w = w
        self.acceleration = acceleration

        self.N_thermal = N_thermal
        
        self.N_sweeps = N_sweeps
        
        self.N_measurement = N_measurement
        
        self.a = a
        
        self.N = N
        
        self.Delta = Delta
        
        self.lattice = np.zeros(N)
        
        self.total_time = N*a
        
        self.sweeps = 0
        
        self.N_tau = N_tau
        self.DH = np.zeros(N_measurement)
        
        self.d_tau = d_tau
        self.accepted = 0
        self.tries = 0
        if N_thermal!= 0:
            print('Termalisation with a = ' + str(self.a) + " N = " +str(self.N)+ " and Accel = " + str(self.acceleration))
            self.thermalize()
            print('----------------------------------')
        self.accepted = 0
        self.tries = 0
        
    def A_k(self,k, M):
        """
        Computes the inverse of the Kernel in Fourier Space
        """
        return self.a/(self.m*(4*np.sin(np.pi*k/self.N)**2+M**2))
    
    def HMC_accel(self):
        """
        Hybrid Monte Carlo with Fourier Acceleration
        """

        N_2 = int(self.N/2)
        k = np.arange(0, self.N) # lattice in Fourier space
        A_k = (self.m/self.a * (4*np.sin(np.pi * k/self.N)**2 + (self.a*self.w)**2) )**(-1)
        self.A = A_k
        PI_std = np.sqrt(self.N /  A_k) 
        Pi = np.random.normal(loc=0, scale=PI_std)

        ps_bar = np.zeros(self.N,dtype=complex)
        ps_bar[0] = Pi[0]
        ps_bar[N_2] = Pi[N_2]


        # unpack remaining PI into Fourier space pi

        ps_bar[1:N_2] = 1/np.sqrt(2) * (Pi[1:N_2] + 1j * Pi[N_2+1:][::-1])

        # impose hermitean symmetry
        ps_bar[N_2+1:] = np.conj(ps_bar[1:N_2][::-1])

        # pi is real by construction
        phi = self.lattice.copy()


        ps = np.real(np.fft.ifft(ps_bar))
        #ps = np.ones(100)
        #phi = np.ones(100)

        H = 1/(2*self.N) *sum(np.abs(np.fft.fft(ps))**2*A_k) +Lattice.action(phi,self.a,self.m,self.w)
        ps_f, phi_f = self.molecular_dynamics_accelerated(phi.copy(),ps.copy(),A_k)


        H_new = 1/(2*self.N) *sum(np.abs(np.fft.fft(-ps_f))**2*A_k) +Lattice.action(phi_f,self.a,self.m,self.w)

        delta_H = H_new- H
        self.delta_H = delta_H

        if delta_H <0 or np.exp(-delta_H) > np.random.uniform(0,1):
            self.accepted += 1
            self.lattice = phi_f.copy()
        self.tries += 1

    def molecular_dynamics_accelerated(self,phi_0,ps_0,A_k):
        """
        Molecular Dynamics with Fourier Acceleration
        """

        p = ps_0 -self.d_tau/2 * (self.m/self.a *(2*phi_0-np.roll(phi_0,1)-np.roll(phi_0,-1))+ self.a*self.m*self.w**2 *phi_0 )
        #print(p)
        p_f = np.fft.fft(p)
        phi =phi_0 +np.real(np.fft.ifft(p_f*A_k))*self.d_tau
        #print(phi)

        for i in range(0,self.N_tau):

            p = p -self.d_tau* (self.m/self.a *(2*phi-np.roll(phi,1)-np.roll(phi,-1))+ self.a*self.m*self.w**2 *phi )
            p_f = np.fft.fft(p)
            phi =phi +np.real(np.fft.ifft(np.multiply(p_f,A_k)) )*self.d_tau

        p = p -self.d_tau/2* (self.m/self.a *(2*phi-np.roll(phi,1)-np.roll(phi,-1))+ self.a*self.m*self.w**2 *phi )
        return p,phi
    

    
    @staticmethod
    def action(lattice,a, m, w):
        """
        Gives the Action of a given lattice.
        """
        
        S =  sum(1/2 *a *((np.roll(lattice,-1)-lattice)/a)**2 +1/2*a*m*w**2*lattice**2)
       
        return S
    
    def HMC(self):
        """
        Hybrid Monte Carlo for the case of no acceleration
        """

        ps = np.array([np.random.normal() for i in range(self.N)])

        H = sum(ps**2)/2 + Lattice.action(self.lattice,self.a,self.m,self.w)

        p_f, x_f = self.molecular_dynamics(ps.copy(),self.lattice.copy())

        H_new = sum(p_f**2)/2 + Lattice.action(x_f,self.a,self.m,self.w)

        self.delta_H = H_new- H

        if self.delta_H <0 or np.exp(-self.delta_H) > np.random.uniform(0,1):
            self.accepted += 1
            self.lattice = x_f.copy()
        self.tries += 1
    

    def molecular_dynamics(self,p_0, x_0):
        """
        Molecular Dynamics for the case of no acceleration
        """
        p = p_0 -self.d_tau/2 * (self.m/self.a *(2*x_0-np.roll(x_0,1)-np.roll(x_0,-1))+ self.a*self.m*self.w**2 *x_0 )
        x = x_0 + self.d_tau*p

        for j in range(1,self.N_tau):

            p = p -self.d_tau * (self.m/self.a *(2*x-np.roll(x,1)-np.roll(x,-1))+ self.a*self.m*self.w**2 *x )
            x = x + self.d_tau*p

        p = p -self.d_tau/2 * (self.m/self.a *(2*x-np.roll(x,1)-np.roll(x,-1))+ self.a*self.m*self.w**2 *x )

        return p, x
        
    
    def thermalize(self):
        """
        Runs the HMC alogrithm for N_thermal times
        """
        with alive_bar(self.N_thermal) as bar:
            for i in range(self.N_thermal):
                if self.acceleration:
                    self.HMC_accel()

                else:
                    self.HMC()
                bar()
    def Calibration_Runs(self, N_runs, N_thermal):
        """
        Runs the HMC for Calibration of the parameters of the algorithms
        """
        with alive_bar(N_runs+N_thermal) as bar:
            for i in range(N_thermal):
                if self.acceleration:
                    self.HMC_accel()
                else:
                    self.HMC()
                bar()
            self.accepted = 0
            self.tries = 0
            for i in range(N_runs):

                if self.acceleration:
                    self.HMC_accel()
                else:
    
                    self.HMC()
                bar()
    

    
    def Get_configurations(self):
        """
        Runs the HMC N_measurment times and records the lattice in each measurement.
        """
    
        results = np.empty((self.N_measurement,self.N))
        with alive_bar(self.N_measurement) as bar:

    
            for i in range(self.N_measurement):
        
                for j in range(self.N_sweeps):
        
                    if self.acceleration:
                        self.HMC_accel()

                    else:
                        self.HMC()
            
            
                results[i] = self.lattice.copy()
                bar()

    
        return(results)
    
    def generate_measurements(self, observable):
        """
        Runs the HMC N_measurment times and records the observable value in each measurement.
        """

        results = [0 for i in range(self.N_measurement)]
        print('Measurements with a = ' + str(self.a) + " N = " +str(self.N)+ " and Accel = " + str(self.acceleration))
        with alive_bar(self.N_measurement) as bar:


            for i in range(self.N_measurement):

                for j in range(self.N_sweeps):

                    if self.acceleration:
                        self.HMC_accel()

                    else:
                        self.HMC()

                if i% 1000 == 0:
                    pass
                self.DH[i]  = self.delta_H 
                results[i] = observable(self.lattice)
                bar()
        print('------------------------------')
        return results

    @staticmethod
    def measure_twopoint(lattice):
        """
        Returns the two point function of a given lattice
        """
        N = len(lattice)

        positions = [i for i in range(N)]
        values = [0 for i in range(positions)]
      
        for i in range(len(positions)):
            values[i] = lattice[0]*lattice[positions[i]]

        return values
    
    @staticmethod    
    def measure_greens(lattice):
        """
        Returns the Greens function of a given lattice
        """

        N = len(lattice)
        
        values = [0 for l in range(N+1)]

        
        for n in range(N+1):
            values[n] = 1/N* sum( np.roll(lattice,-n)*lattice)
        

        return values

    @staticmethod
    def measure_position(lattice):
        """
        Returns the average position of a given lattice
        """
        return np.average(lattice)
    @staticmethod
    def measure_sq_position(lattice):
        """
        Returns the average sq position of a given lattice
        """
        return np.average([x**2 for x in lattice])
    
    @staticmethod
    def save_measurements(N, a, N_thermal, N_meausre, N_sweeps, measurements, observable_name, file_name):
        """
        Option to save files
        """
        file = open(file_name,'a')
        
        file.write("########################"+ '\n')
        
        file.write("Description:" + '\n')

        file.write('The Data file that comes with this one stores the measremts of the observable ' + observable_name+ ". The following lines give a breif desription of the simulation used to draw them.")
        
        file.write('The lattice has '+ str(N)+ "sites and a is chosen to be " + str(a) + ' in this simulation'+ '\n')
        
        file.write(str(N_thermal) + ' sweeps have been at the beginning to thermalize the lattice'+ '\n')
        
        file.write(str(N_meausre) + ' Configurations of the lattice have been saved'+ '\n')
        
        file.write(str(N_sweeps)+ ' sweeps have been performed between each saved configuration to minimize the autocorrelation'+ '\n')
        
        file.write("########################"+ '\n')
        
        file.close()

        np.save(file_name,measurements)
