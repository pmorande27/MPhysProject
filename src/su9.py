



import numpy as np
import cProfile
import time
from Chiral import Chiral
import matrix_routines as Mat
import plotting
import Chiral_run
import matplotlib.pyplot as plt
import processing
import Greens

def main():

    SU = 9
    order = 10
    N_order = 10
    betas2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.6,1.7,1.8,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,8.0]
    N = [16 for i in range(len(betas2))]
    N_thermal = [10**4 for i in range(len(N))]
    N_tau,_ =10,1/10
    accel = True
    N_measure = [10**4 for i in range(len(N))]
    N_tau = 30
    for i,beta in enumerate(betas2):
        #N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau,True,accel)
        
        #Chiral_run.measure(beta,N[0],SU,order,N_order,N_measure[0],N_thermal[0],lambda U: -Chiral.action(U,beta,2)/(2*SU*N[0]**2*beta),"Action",True,accel)
        processing.process_Action(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel) 
        pass
    model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    #plotting.plot_generic_data_beta(betas2,model_params,'Action','e',accel=accel)
    plotting.plot_e_desinty(betas2,model_params,accel=accel)
    
    #iats, beta = processing.sus_iat_process(betas1,16,SU,order,N_order,10000,1000)
    #plt.plot(beta,iats)
    #plt.show()
main()


