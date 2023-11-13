import numpy as np
import Stats
def Greens_mom(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.array([[Stats.Stats(data[i][j]).estimate()[0] for j in range(N)]for i in range(N)])
    Greens_zero_mom = np.zeros((N+1))
    for i in range(N):
        Greens_zero_mom[i] = np.sum(values[i])
    Greens_zero_mom[N] = Greens_zero_mom[0]
    observable_name = 'Greens 0 Mom'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    np.save(file_name,Greens_zero_mom/Greens_zero_mom[0])

def process_ww_corr(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'ww corr'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i]= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'

    np.save(file_name,(values,values_err))