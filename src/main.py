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

    N = [30]
    SU = 3
    betas1 = [1.5]
    order = 10
    N_order = 10
    N_tau = 60
    N_thermal = 10**4
    N_measure = 10**4
    for i,beta in enumerate(betas1):
        #N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau)
        Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure,N_thermal,lambda U:Chiral.Measure_G(U,SU),"Greens")
        pass
    model_params = {'n_length': N[0],'su_parameter': SU, 'order': order, 'n_measure': N_measure, 'n_thermal': N_thermal, 'n_order':N_order }
    #[processing.process_Action(betas1[i],N[0],SU,order,N_order,N_measure,N_thermal) for i in range(len(betas1))]
    #plotting.plot_e_desinty(betas=betas1,model_params=model_params)
    processing.process_Greens(betas1[0],N[0],SU,order,N_order,N_measure,N_thermal)
    #Greens.Greens_mom(betas1[0],N[0],SU,order,N_order,N_measure,N_thermal)
    #plotting.plot_Greens_0_mom(betas1[0],N[0],SU,order,N_order,N_measure,N_thermal)

  
main()

