import numpy as np
import cProfile
import time
from Chiral import Chiral
import matrix_routines as Mat
import plotting
import Chiral_run
import Stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import processing

def main():

    N = [64]
    SU = 2
    betas1 = [0.8667]
    order = 0
    N_order = 0
    N_tau = 60
    N_thermal = 10**4
    N_measure = 10**6
    for i,beta in enumerate(betas1):
        #N_tau = Chiral_run.calibration(beta,N[i],SU,order,N_order,N_tau)
        #Chiral_run.measure_func_1D(beta,N[i],SU,order,N_order,N_measure,N_thermal,lambda U:Chiral.Measure_ww_corr(U,SU),"ww corr")
        pass
    model_params = {'n_length': N,'su_parameter': SU, 'order': order, 'n_measure': N_measure, 'n_thermal': N_thermal, 'n_order':N_order }
    #print(sus2(1.08,N,SU,order,N_order,N_measure,N_thermal))
    #observable_name = 'Susceptibility'
    #file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(1.08) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    #print(Stats.Stats(np.load(file_name)).estimate()[0])
    #print(corr_len(1.08,N,SU,order,N_order,N_measure,N_thermal))
    #Greens_mom(0.8667,N,SU,order,N_order,N_measure,N_thermal)
    processing.process_ww_corr(betas1[0],N[0],SU,order,N_order,N_measure,N_thermal)
    plotting.Greens_mom_2(betas1[0],N[0],SU,order,N_order,N_measure,N_thermal)
    
    """observable_name = 'Greens 0 Mom'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(0.8667) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    
    ydata = np.load(file_name)
    xdata = np.arange(N+1)
    def model(t, dE):
        L = N
        return  (np.cosh(1/dE*(t-L/2))-1 ) / (np.cosh(1/dE*L/2) -1)
    a =curve_fit(model,xdata,ydata,absolute_sigma=True,)
    xs = [model(d,a[0][0]) for d in xdata]
    print(a[0][0])
    plt.plot(xdata,xs)

    plt.plot(xdata,np.load(file_name),label='holaaa')
    plt.legend()
    plt.show()"""
    

    #plotting.plot_e_desinty([0.1,0.2,0.3], model_params)
    #plotting.plot_generic_data_beta(betas2, model_params,'Susceptibility','$\chi$')
    #print(load_calibration(7.0,N,SU,order,N_order))
    #a = Chiral_run.estimate_2D_func(1.08,N,SU,order,N_order,N_measure,N_thermal,lambda U:Chiral.Measure_G(U,SU),"Greens")
    #b = np.fft.fft2(a)
    #print(np.sqrt(1/(4*np.sin(np.pi/16)**2  ) * (b[0,0]/b[0,1]-1)).real)
def sus3(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 2'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    print(data.shape)
    values = np.array([[[[Stats.Stats(data[i][j][m][n]).estimate()[0] for n in range(N)]for m in range(N)]for j in range(N)] for i in range(N)])
    print('a')
    result = 0
    sus = 0
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    #result += (i**2+j**2)*values[i][j]
                    sus += values[i][j][m][n]

    return sus/N**2
def corr_len_2(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 2'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    print(data.shape)
    values = np.array([[[[Stats.Stats(data[i][j][m][n]).estimate()[0] for n in range(N)]for m in range(N)]for j in range(N)] for i in range(N)])
    result = 0
    sus = 0
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    result += ((i-m)**2+(j-n)**2)*values[i][j][m][n]
                    sus += values[i][j][m][n]
    result = result/N**2
    sus = sus /N**2
    result = result/(4*sus)
    return (result**0.5,sus)
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

    
    
def sus2(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    print(data.shape)
    values = np.array([[Stats.Stats(data[i][j]).estimate()[0] for j in range(N)]for i in range(N)])
    result = 0
    sus = 0
    for i in range(N):
        for j in range(N):
            result += (i**2+j**2)*values[i][j]
            sus += values[i][j]

    return (result/(4*sus))**0.5,sus
def corr_len(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    print(data.shape)
    values = np.array([[Stats.Stats(data[i][j]).estimate()[0] for j in range(N)]for i in range(N)])
    result = 0
    sus = 0
    Gf = np.fft.fft2(values)
    p = 2*np.pi/N
    """Gf1 =np.sum(np.array([ [values[i,j] * np.exp(1j*(p)*j) for j in range(N)]for i in range(N)]))
    #print(np.sum(Gf1))
    Gf0 = np.sum(np.array([[values[i,j] for j in range(N)] for i in range(N)]))
    """
    return (1/(4*np.sin(np.pi/N)*np.sin(np.pi/N))*(Gf[0,0]/Gf[0,1]-1))**0.5,Gf[0,0]
main()

