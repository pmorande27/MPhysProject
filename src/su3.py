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
from scipy.optimize import curve_fit

def main():
    betas2= [0.1,0.2,0.3,0.4,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    N = [16 for b in betas2]
    SU = 3
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
        #processing.process_CV_jack(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
    model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    #plotting.plot_e_desinty(betas2,model_params,accel=accel)
    #print(2* N[0]**2*SU)
    N_orders = [0,10]
    orders = [0,10]
    N_measures =[10**6,10**5]
    N_thermals = [10**4,10**4]
    plotting.plot_CV_mult([betas2,betas2], model_params,[2,3],N_measures,N_thermals, orders,N_orders,accel = True)     
main()


