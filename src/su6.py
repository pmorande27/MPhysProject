



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

    """SU = 6
    order = 10
    N_order = 10
    betas2 = [3.7199999999999998 ]
    N = [54]
    N_thermal = [10**4 for i in range(len(N))]
    N_tau,_ =10,1/10
    accel = True
    N_measure = [10**5 for i in range(len(N))]
    
    for i,beta in enumerate(betas2):
        print(beta)
        #N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau,True,accel)
        #Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr",True,accel)
        Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr",True,accel)
        #processing.process_ww_corr(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel) 
        #processing.process_ww_corr_state_two(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel) 
        #Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_G(U,SU),"Greens",True,accel)
        #processing.process_Greens(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U: -Chiral.action(U,beta,2)/(2*SU*N[0]**2*beta),"Action",True,accel)
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
    #model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    #plotting.plot_e_desinty(betas2,model_params,accel=accel)
    
    #iats, beta = processing.sus_iat_process(betas1,16,SU,order,N_order,10000,1000)
    #plt.plot(beta,iats)
    #plt.show()"""
    betas2= [2.3]
    N = [16 for b in betas2]
    SU = 6
    Ns= [30,42,48,90,120]
    accel = True
    order = 10
    #processing.process_mass_cosh_jack(beta,N,SU,order,N_order,N_measure,N_thermal,accel=accel)
    #plotting.plot_iats_susceptibility(beta,N,SU,order,N_order,10**5,10**4)
    
    N_order = 10
    N_measure = [10**5 for b in betas2]
    N_thermal = [10**4 for b in betas2]
    #plotting.plot_cosh_mass_jack_def(1.2667,400,SU,order,N_order,10**5,10**4,range(15,20),[30],accel=accel)
    #lotting.correlation_length_fit_two_params(1.2667,400,SU,order,N_order,10**5,10**4,30,400,accel=accel)
    #plotting.plot_SU2_massoverlambda(betas,Ns,order,N_order,N_measure,N_thermal)
    #processing.process_mass_cosh_jack(beta,N,SU,order,N_order,N_measure,N_thermal,accel=accel)
    #plotting.plot_iats_susceptibility(betas,Ns,3,order,N_order,10**5,10**4)
    #plotting.plot_cosh_mass_jack_def(beta,N,SU,order,N_order,10**5,10**4,[20,21,22],[40,41,42],accel=accel)
    for i,beta in enumerate(betas2):
        #N_tau = Chiral_run.calibration(beta,N[0],SU,order,N_order,N_tau,True,accel)
        Observable_1 = lambda U: (Chiral.action(U,beta,2)/beta)**2
        Observable_2 = lambda U: (Chiral.action(U,beta,2)/beta)
        Chiral_run.measure_two_obs(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],Observable_1,Observable_2,"CV_data",True,accel)
        #processing.process_CV_jack(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
    
main()


