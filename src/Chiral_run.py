from Chiral import Chiral
import numpy as np
from Stats import Stats
import matplotlib.pyplot as plt
def strong_coupling(beta):
    return 1/2*beta+1/6*beta**3+1/6*beta**6
def weak_coupling(beta):
    Q_1 = 0.0958876
    Q_2 = -0.067
    return 1- 3/(8*beta)*(1+1/(16*beta)+(1/64+3/16*Q_1+1/8*Q_2)*1/beta**2)
def calibration(beta, N, SU, N_tau_guess = 2):
    N_tau = N_tau_guess
    print('Calibration with beta = ' + str(beta) + " N = " +str(N)+ " SU = " + str(SU) )
    up = 0.9
    low = 0.75
    max_count = 10
    results = [0 for i in range(max_count)]
    for i in range(max_count):
        epsilon = 1/N_tau
        calibration_runs = 10**3
        lat = Chiral(N, beta, 0,0,1,epsilon, N_tau, SU)
        lat.Calibration_Runs(calibration_runs, 1000)
        rate = lat.accepted/lat.tries
        d_rate = 0.65-rate
        results[i] = (rate-up,N_tau)
        print(rate,N_tau)
        

        new_N = int(np.rint(N_tau*(1+d_rate)))
        if rate <=up and rate >= low:
            file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)
            np.save(file_name, [N_tau,1/N_tau])
            print("-----------------")
            print(rate, N_tau)
            return N_tau
        if new_N == N_tau:
            if d_rate <0:
                new_N -= 1
            else:
                new_N +=1
        
        N_tau = new_N
       

    print("-----------------")
    print("Calibration Unsucessful, better run:")
    results_abs = [(abs(x),y) for (x,y) in results]
    d_rate, N_tau = min(results_abs)
    d_rate_2 = lookup(d_rate,N_tau,results)
    rate = (d_rate_2+up)*100
    print(rate,N_tau)
    file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)
    np.save(file_name, [N_tau,1/N_tau])
    return N_tau
def load_calibration(beta, N, SU):
    file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+'.npy'
    values = np.load(    file_name)
    return int(values[0]),values[1]
def lookup(d_rate,N_tau,results):
    for (x,y) in results:
        if abs(x) == d_rate and y == N_tau:
            return x
def measure_action(beta,N,SU):
    file_name = "ChiralResults/Action 2 beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+'.npy'
    N_tau, epsilon = load_calibration(beta,N,SU)
    model = Chiral(N,beta,10**3,1000,1,epsilon,N_tau,2,1,order=10)
        
    results = model.generate_measurements(lambda U: -Chiral.action(U,beta)/(4*N**2*beta))
    np.save(file_name,results)


def plot_e_desinty(betas,N, SU):
    result = np.zeros(len(betas))

    error = np.zeros(len(betas))
    result_2 = np.zeros(len(betas))

    error_2 = np.zeros(len(betas))

    

    for i in range(len(betas)):

        file_name = "ChiralResults/Action"   + " beta = " + str(betas[i]) +" N = "+str(N)+ " SU = " + str(SU)+".npy"

        file_name_2 = "ChiralResults/Action 2"   + " beta = " + str(betas[i]) +" N = "+str(N)+ " SU = " + str(SU)+".npy"

        values = np.load(file_name)
        values_2  = np.load(file_name_2)
        result[i],error[i] = Stats(values).estimate()
        result_2[i],error_2[i] = Stats(values_2).estimate()

        

    ax = plt.subplot(111)
    
    ax.errorbar(x=betas,y=result,yerr= error,fmt=".k",label = 'Data')
    ax.errorbar(x=betas,y=result_2,yerr= error_2,fmt=".k",label = 'Data with Taylor',color = 'g')

    betas_list_1 = np.linspace(0,1,100)
    strong = [strong_coupling(beta) for beta in betas_list_1]
    plt.plot(betas_list_1,strong,'b')
    betas_list_2 = np.linspace(0.5,3.0,100)
    weak = [weak_coupling(beta) for beta in betas_list_2]
    plt.plot(betas_list_2,weak,'r')




    plt.xlim(0,3.1)
    
    plt.legend()
    
    plt.xlabel('$\\beta$')
    
    plt.ylabel('$e$')
    
    ax.spines[['right', 'top']].set_visible(False)
    

    plt.show()
def main():
    N = 16
    SU = 2
    betas1 = [2.9,3.0]
    betas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    N_tau = 18
    for beta in betas1:
        #N_tau = calibration(beta,N,SU,N_tau)
        #measure_action(beta,N,SU)
        pass
    plot_e_desinty(betas,N,SU)

main()
