import numpy as np
import matrix_routines as Mat
def main():
    a = np.load('beta_1_0.npy')

    values = np.concatenate((a[1],a[1][::-1][1:]))

    errs = np.concatenate((a[2],a[2][::-1][1:]))
    beta = 1.0
    
    N = 96
    SU = 2
    order = 0
    N_order = 0
    N_measure = 10**5
    N_thermal = 10**3
    observable_name = 'ww corr'
    print(values)
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(values,errs))
main()
