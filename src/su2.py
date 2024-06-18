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

    SU = 2
    betas2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    #betas2 = [0.1,0.2]
    #N = [18,24,30,36,42,48,81,90,120]
    #betas1 = [0.6,0.6667,0.7333,0.8,0.9333]
    #cor_len =[ 1.425,1.776,2.268,2.958, 3.987,5.542,7.839]
    N = [16 for b in betas2]
    betas1 = [0.8]
    order = 0
    N_order = 0
    #N_tau,e=Chiral_run.load_calibration(betas1[0],N[0],SU,10,10)
    N_thermal = [10**4 for i in range(len(N))]
    N_tau,_ =1,1
    accel = True
    N_measure = [10**4 for i in range(len(N))]
    #plotting.plot_iats(betas1,N,SU,order,N_order,N_measure[0],N_thermal[0],corr_lens=cor_len)
    for i,beta in enumerate(betas2):
        N_tau = Chiral_run.calibration(beta,N[0],SU,order,N_order,N_tau,True,accel)
        #Observable_1 = lambda U: (Chiral.action(U,beta,2)/beta)**2
        #Observable_2 = lambda U: (Chiral.action(U,beta,2)/beta)
        #Chiral_run.measure_two_obs(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],Observable_1,Observable_2,"CV_data",True,accel)
        #processing.process_CV_jack(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U: (-Chiral.action(U,beta,2)/beta)**2/(2*SU*N[0]**2),"Action Sq",True,accel)
        #Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U: (-Chiral.action(U,beta,2)/beta)/(2*SU*N[0]**2),"Action",True,accel)
        #processing.process_Action(betas2[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_Action_sq(betas2[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_CV(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
       
        #Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_G(U,SU),"Greens",True,accel)
        #processing.process_Greens(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_Action(betas2[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #print(Greens.second_moment_correletion_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel))
        #Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr",True,accel)
        #processing.process_ww_corr(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_Generic_observable(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],"Susceptibility",accel)
        #processing.process_iats(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #Chiral_run.measure(beta,N[0],SU,order,N_order,N_measure[0],N_thermal[0],lambda U: -Chiral.action(U,beta,2)/(2*SU*N[0]**2*beta),"Action",True,accel)
        #processing.process_ww_corr(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i])
        #plotting.Greens_mom_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i])
        #processing.reprocess_ww(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],5*N_thermal[i])
        #plotting.correlation_length_manual_fit(betas1[i],N[i],SU,order,N_order,N_measure[i]-5*N_thermal[i],6*N_thermal[i])
        #plotting.correlation_length_manual_fit(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i])
        #=print(Greens.second_moment_correletion_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel))
        #plotting.correlation_length_manual_fit(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],lower=2,upper = 25,accel=accel)
        #plotting.correlation_length_automatic_fit(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],lower=2,accel=accel)
        #plotting.correlation_length_automatic_fit(betas1[i],N[i],SU,order,N_order,N_measure[i]-5*N_thermal[i],6*N_thermal[i])
        #Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_G(U,SU),"Greens")
        #processing.process_Greens(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i])
        #plotting.mass_plot(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i])
        #print(Greens.second_moment_correletion_length_two(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i]))

        #print(Greens.second_moment_correletion_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i]))
        #processing.process_Action(betas1[i],N[0],SU,order,N_order,N_measure[0],N_thermal[0],accel=accel)
        pass
    model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    plotting.plot_e_desinty(betas2,model_params,accel=accel)
    #print(2* N[0]**2*SU)
    
    
    plotting.plot_CV(betas2, model_params, accel = True)
    
    #iats, beta = processing.sus_iat_process(betas1,16,SU,order,N_order,10000,1000)
    #plt.plot(beta,iats)
    #plt.show()
main()


