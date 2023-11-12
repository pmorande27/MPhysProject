import numpy as np
import cProfile
import time
from Chiral import Chiral
import matrix_routines as Mat
import plotting
import Chiral_run
import Stats

def main():

    N = 48
    SU = 3
    betas1 = [1.62]
    betas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0]
    order = 10
    betas2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
    N_order = 10
    N_tau = 20
    N_thermal = 10**3
    N_measure = 1000
    for beta in betas1:
        N_tau = Chiral_run.calibration(beta,N,SU,order,N_order,N_tau)
        Chiral_run.measure_func_4D(beta,N,SU,order,N_order,N_measure,N_thermal,lambda U:Chiral.Measure_G(U,SU),"Greens")
        pass
    model_params = {'n_length': N,'su_parameter': SU, 'order': order, 'n_measure': N_measure, 'n_thermal': N_thermal, 'n_order':N_order }
    #print(sus2(1.08,N,SU,order,N_order,N_measure,N_thermal))
    observable_name = 'Susceptibility'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(1.08) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    #print(Stats.Stats(np.load(file_name)).estimate()[0])
    print(corr_len(1.62,N,SU,order,N_order,N_measure,N_thermal))

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
def main2():
    SU = 3
    N = 2
    Identity = np.random.uniform(size=(N,N,SU,SU))+1j*np.random.uniform(size=(N,N,SU,SU))
   
   
    print(np.einsum('ijkl,mnlk->',Identity,Mat.dagger(Identity)).real/(3*N**2))
    print(susceptibility_naive(Identity))
    print(1/SU *np.einsum('ijkl,lk->ij',Identity,Mat.dagger(Identity)[0,0]).real)

def susceptibility_naive(phi):
    N = len(phi)
    G = np.zeros((N,N),complex)
    for i in range(N):
        for j in range(N):
            phi_y = np.roll(phi, shift=(i,j), axis=(0,1))
            A = Mat.multiply_matrices(phi, Mat.dagger(phi_y))
            G += np.einsum('ijkk->ij',A+Mat.dagger(A))

    Chi = np.sum(G) / (2*3*N**2)

    return Chi
#main2()