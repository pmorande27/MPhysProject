import numpy as np
import cProfile
import time
from Chiral import Chiral
import matrix_routines as Mat
import plotting
import Chiral_run

def main():

    N = 16
    SU = 3
    betas1 = [9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0]
    betas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,11.0,12.0,13.0]
    order = 10
    betas2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
    N_order = 10
    N_tau = 79
    N_thermal = 10**3
    N_measure = 10**4
    for beta in betas2:
       #N_tau = Chiral_run.calibration(beta,N,SU,order,N_order,N_tau)
        #Chiral_run.measure(beta,N,SU,order,N_order,N_measure,N_thermal,lambda U: -Chiral.action(U,beta,2)/(2*SU*N**2*beta),"Action")
        pass
    model_params = {'n_length': N,'su_parameter': SU, 'order': order, 'n_measure': N_measure, 'n_thermal': N_thermal, 'n_order':N_order }
    plotting.plot_e_desinty(betas2, model_params)
    #plotting.plot_generic_data_beta(betas2, model_params,'Action','e')
    #print(load_calibration(7.0,N,SU,order,N_order))


main()
def main2():
    SU = 4
    N = 4
    Identity = np.zeros((N,N,SU,SU))
    for i in range(N):
        for j in range(N):
            Identity[i,j] = np.eye(SU)
    D =np.random.uniform(size = (N,N,SU,SU))+ 1j*np.random.uniform(size = (N,N,SU,SU))
    E = Mat.reunitarisation(D.copy(),4)
    np.testing.assert_allclose(Mat.determinant(E),np.ones((N,N)))
    np.testing.assert_allclose(Mat.multiply_matrices(E,Mat.dagger(E)),Identity+np.zeros((N,N,SU,SU),complex),atol=10**-14,rtol=10**-16)
#main2()