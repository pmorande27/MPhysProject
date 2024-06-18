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
import Stats

def main():

    """SU = 4
    #N = [16]
    #betas1 = [round(i,1) for i in np.arange(0.1,9.8,0.1)]
    #N = [16 for b in betas1]
    #cor_len =[ 1.425,1.776,2.268,2.958, 3.987,5.542,7.839]
    #betas1 = [1.2]
    #N = [16 for i in betas1]
    def rounding(x):
        y = str(x)
        l = ""
        flag = False
        for (i,val) in enumerate(y):
            
            
            if i != len(y)-1:
                if val == "0" and y[i+1] == "0":
                    flag = True
                    break
            l += val
        
        if not flag:
            l = ""
            for (i,val) in enumerate(y):
                
                if i != len(y)-1:
                    if y[i+1] == "9" and y[i+2] == "9":
                        l += str(int(val)+1)
                        break
                l += val
        return float(l)
    
    betas1 = [ rounding(8*i) for i in [0.18,0.225,0.25,0.27,0.27,0.27,0.29,0.30,0.315]]
    mass =1/np.array( [1.003,1.87,3.027,4.79,4.79,4.81,7.99,10.41,15.5])
    N = [18,24,30,36,42,48,82,90,120]
    print(betas1)
    order = 10
    N_order = 10
    #N_tau,e=Chiral_run.load_calibration(betas1[0],N[0],SU,10,10)
    N_thermal = [10**4 for i in range(len(N))]
    accel = True
    N_measure = [10**5 for i in range(len(N))]

    #N = 16
    #beta = 1.2
    #N_taus = [5,6,7,8,9,10,11,12,13]
    N_tau =1

    #Chiral_run.exponential_measurements(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],N_taus,accel=accel)
    #plotting.plot_accH(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],accel=accel)
    #plotting.plot_exponential(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],accel=accel)
    for i,beta in enumerate(betas1):
        N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,1,True,accel)
        Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr",True,accel)
        Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr_state_two(U,SU),"state 2",True,accel)
        Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_G(U,SU),"Greens",True,accel)
        pass"""
    betas2= [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0]
    N = [16 for b in betas2]
    SU = 4
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
        #Chiral_run.measure_two_obs(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],Observable_1,Observable_2,"CV_data",True,accel)
        processing.process_CV_jack(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
    model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    #plotting.plot_e_desinty(betas2,model_params,accel=accel)
    #print(2* N[0]**2*SU)
    
    
    plotting.plot_CV(betas2, model_params, accel = True)    
  
main()


