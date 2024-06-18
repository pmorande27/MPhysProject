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

    SU = 3
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
    
    betas1 = [    1.08  ]
    N = [18]
    order = 10
    N_order = 10
    #N_tau,e=Chiral_run.load_calibration(betas1[0],N[0],SU,10,10)
    N_thermal = [10**4 for i in range(len(N))]
    accel = True
    N_measure = [10**5 for i in range(len(N))]
    mass = np.zeros(len(betas1))
    for i,beta in enumerate(betas1):
        observable_name = "mass state 1"
        file_name ="ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "   + str(N_measure[i])+" N Thermal = "  + str(N_thermal[i])+' Accel.npy'
            
        x,err_x = np.load(file_name)
        mass[i] = x
    
    #tminss = [[2,3,4],[2,3,4],[3,4,5],[5,6,7,8]]
    #tmaxss = [[7,8,9],[10,11],[14,15,16,17],[18,19,20]]

    #tminss = [[1,2,3]]
    #tmaxss = [[5]]
    #tminsse = [[1,2,3]]

    #tmaxss_cor = [[1] ]
    #tminss_cor =[[0]]

    #N = 16
    #beta = 1.2
    #N_taus = [5,6,7,8,9,10,11,12,13]

    #Chiral_run.exponential_measurements(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],N_taus,accel=accel)
    #plotting.plot_accH(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],accel=accel)
    #plotting.plot_exponential(beta,N,SU,order,N_order,N_measure[0],N_thermal[0],accel=accel)
    for i,beta in enumerate(betas1):
        Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_Sus(U,SU),"Susceptibility",True,True,mass[i])
        
        #Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_Sus(U,SU),"Susceptibility",True,False)
        #processing.process_iats(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=True)
        #processing.process_iats(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=False)
        #Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr",True,accel)
        #processing.process_mass_log_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_mass_cosh_jack_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_mass_cosh_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_second_moment_correlation_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_ww_corr_state_2_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_ww_corr_jack_s(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_mass_log_jack_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #plotting.correlation_length_fit_two_params_2_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[15,16],[25,26,27,28],accel=accel)
        #plotting.plot_cosh_mass_jack_def(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[10,11,12,13,14,15],[23,24,25,26],accel=accel)
        #plotting.plotlog_log_mass_jack_def_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[20,21,22],[28,27],accel=accel)

        #plotting.plot_cosh_mass_jack_def_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[16,17],[22],accel=accel)
        """observable_name = "mass state 1"
        file_name ="ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "   + str(N_measure[i])+" N Thermal = "  + str(N_thermal[i])+' Accel.npy'
        
        x,err_x = np.load(file_name)
        #x,err_x = 1/x, err_x/x**2 
        observable_name = "mass state 2"
        file_name ="ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "   + str(N_measure[i])+" N Thermal = "  + str(N_thermal[i])+' Accel.npy'
        y,err_y = np.load(file_name)
        #print(y)
        z = y/x
        print(y)
        err_z = np.sqrt((err_y/x)**2+(err_x*y/x**2)**2)

        print(z,err_z )"""

        #plotting.plot_accH(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)

        #plotting.plotlog_log_mass_jack_def_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[5],[7,8,9],accel=accel)
        #plotting.plot_cosh_mass_jack_def_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[11,12,13],[19,20,21],accel=accel)
        #plotting.plot``
        #processing.process_mass_cosh_jack_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        
       # plotting.plot_cosh_mass_jack_state_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],11,7,accel=accel)
        #processing.process_ww_corr_state_2_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #plotting.correlation_length_fit_two_params_2_state_2_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],6,14,accel=accel)
        #plotting.correlation_length_fit_two_params_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[1],[7],accel=accel)
        #plotting.plot_log_mass_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],4,2,accel=accel)
        #plotting.plot_log_mass_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],4,2,accel=accel)
        #plotting.plot_log_mass_jack_def(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[8,9,10,11,12],[40],accel=accel)
        #print('-------------------------')

        #plotting.plot_cosh_mass_jack_def(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[15,16,17,18,19,20,21,22,23,24],[80],accel=accel)
        #print('-------------------------')
        #plotting.correlation_length_fit_two_params_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],range(2,50),[100,101,102,103],accel=accel)
        #print(np.load("ChiralResults/Processed/second moment correlation length/second moment correlation length beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure[i])+" N Thermal = "  + str(N_thermal[i])+" Accel.npy"))
        #processing.process_ww_corr_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_second_moment_correlation_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #print(Greens.second_moment_correletion_length(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel))
        #N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau,True,accel)
        #Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_ww_corr_state_two(U,SU),"state 2",True,accel)
        #Chiral_run.measure_func_2D(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U:Chiral.Measure_G(U,SU),"Greens",True,accel)
        #processing.process_Greens(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_ww_corr_state_two(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #plotting.correlation_length_def_fit_state_two(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],tminss_cor[i],tmaxss_cor[i],accel=accel)
        #N_tau  = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau,True,accel)
        #Chiral_run.measure(beta,N[i],SU,order,N_order,N_measure[i],N_thermal[i],lambda U: -Chiral.action(U,beta,2)/(2*SU*N[i]**2*beta),"Action",True,accel)
        #processing.process_ww_corr(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #processing.process_Action(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],accel=accel)
        #plotting.cor_plot_def_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],tminss[i],tmaxss[i],accel=accel)
        #plotting.mass_plot_def(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],tminss[i],tmaxss[i],tminss[i],accel=accel)
        #plotting.mass_plot_def_jack(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],tminss[i],tmaxss[i],tminss[i],accel=accel)
        #plotting.mass_plot_def_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],1,10,accel=accel)
        #plotting.correlation_length_def_fit(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],tminss_cor[i],tmaxss_cor[i],accel=accel)
        #plotting.correlation_length_fit_two_params_2(betas1[i],N[i],SU,order,N_order,N_measure[i],N_thermal[i],[1,2,3,4],[7,8,9],accel=accel)
        pass
    #model_params = {'n_length':N[0],'su_parameter':SU,'order':order,'n_order':N_order,'n_measure':N_measure[0],'n_thermal':N_thermal[0]}
    #plotting.plot_e_desinty(betas1,model_params,accel=accel)
    
main()


