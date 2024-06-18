import numpy as np
import plotting
def load_corr_len_function(file_name,beta,N,mode):
    a = np.load(file_name)
    if mode == 1:
        values = np.concatenate((a[1],a[1][::-1][1:]))

        errs = np.concatenate((a[2],a[2][::-1][1:]))
    else:
        values = np.concatenate((a[1],[0],a[1][::-1][1:],[a[1][0]]))

        errs = np.concatenate((a[2],[0],a[2][::-1][1:],[a[2][0]]))
    SU = 2
    order = 0
    print(len(values))
    N_order = 0
    N_measure = 10**5
    N_thermal = 2*10**3
    observable_name = 'ww corr'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(values,errs))
def plot_julian_corr_len(beta,N):
    SU = 2
    order = 0
    N_order = 0
    N_measure = 10**5
    N_thermal = 2*10**3
    observable_name = 'ww corr'
    plotting.correlation_length_def_fit_Julian(beta,N,SU,order,N_order,N_measure,N_thermal,[0],[140],accel=True)
    plotting.mass_plot_def(beta,N,SU,order,N_order,N_measure,N_thermal,[10],[140],[0],accel=True)
    plotting.correlation_length_def_fit(beta,N,SU,order,N_order,N_measure,N_thermal,[10,11,12,13,14,15],[140],accel=True)
    
    


file_name = 'beta_1_2667.npy'
load_corr_len_function(file_name,1.2667,400,1)
plot_julian_corr_len(1.2667,400)

