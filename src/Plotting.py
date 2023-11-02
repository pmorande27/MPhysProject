import numpy as np
from Stats import Stats
import matplotlib.pyplot as plt

def plot_e_desinty(betas,N, SU,order, N_order, N_measure, N_thermal):
    result = np.zeros(len(betas))

    error = np.zeros(len(betas))
    result_2 = np.zeros(len(betas))

    error_2 = np.zeros(len(betas))

    

    for i in range(len(betas)):
        if SU ==2:
            file_name = "ChiralResults/Action/Action beta = " + str(betas[i]) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(0)+" N Order = "  + str(0)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            values = np.load(file_name)
            result[i],error[i] = Stats(values).estimate()


        file_name_2 =  "ChiralResults/Action/Action beta = " + str(betas[i]) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'


        values_2  = np.load(file_name_2)
        result_2[i],error_2[i] = Stats(values_2).estimate()

        

    ax = plt.subplot(111)
    if SU == 2:
        ax.errorbar(x=betas,y=result,yerr= error,fmt="xk",label = 'Data With Exact Exponetnial')
        betas_list_1 = np.linspace(0,1,100)
        strong_2s = [strong_coupling_alt(beta) for beta in betas_list_1]

        strong = [strong_coupling(beta) for beta in betas_list_1]
        plt.plot(betas_list_1,strong_2s,'b',label='Strong Coupling Alt')

        plt.plot(betas_list_1,strong,'b',label='Strong Coupling')
        betas_list_2 = np.linspace(0.5,4.0,1000)
        weak = [weak_coupling(beta) for beta in betas_list_2]
        plt.plot(betas_list_2,weak,'r',label= 'Weak Coupling')
    
    ax.errorbar(x=betas,y=result_2,yerr= error_2,fmt="xg",label = 'Data with Taylor of Order = ' + str(order)+ " and N_order = 10")
    if SU == 3:
        betas = np.linspace(0,2,100)
        values = [strong_2(beta) for beta in betas]
        plt.plot(betas,values,'b',label = 'Strong Coupling Expansion')
        betas = np.linspace(1,4,100)
        values_w= [weak_2(beta) for beta in betas]
        plt.plot(betas,values_w,'r',label = 'Weak Coupling Expansion')
    if SU== 4:
        betas = np.linspace(1,8,100)
        values_w= [weak_3(beta) for beta in betas]
        plt.plot(betas,values_w,'r',label = 'Weak Coupling Expansion')

        betas = np.linspace(0,3,100)
        strong= [strong_3(beta) for beta in betas]
        plt.plot(betas,strong,'b',label = 'Strong Coupling Expansion')

    



    plt.xlim(0,7.1)
    
    plt.legend()
    
    plt.xlabel('$\\beta$')
    
    plt.ylabel('$e$')
    
    ax.spines[['right', 'top']].set_visible(False)
    
    plt.savefig('ChiralResults/Plots/Energy_density_SU_' +str(SU)+"_N_"+str(N)+'Order_'+str(order)+'N_order_'+str(N_order)+'.svg')
    plt.show()
def plot_sus(betas,N, SU):
    result = np.zeros(len(betas))

    error = np.zeros(len(betas))
    result_2 = np.zeros(len(betas))

    error_2 = np.zeros(len(betas))

    

    for i in range(len(betas)):

        #file_name = "ChiralResults/Action/Action"   + " beta = " + str(betas[i]) +" N = "+str(N)+ " SU = " + str(SU)+".npy"

        file_name_2 = "ChiralResults/Susceptibility/Susceptibility 2"   + " beta = " + str(betas[i]) +" N = "+str(N)+ " SU = " + str(SU)+".npy"

        #values = np.load(file_name)
        values_2  = np.load(file_name_2)
        #result[i],error[i] = Stats(values).estimate()
        result_2[i],error_2[i] = Stats(values_2).estimate()

        

    ax = plt.subplot(111)
    
    #ax.errorbar(x=betas,y=result,yerr= error,fmt="xk",label = 'Data')
    ax.errorbar(x=betas,y=result_2,yerr= error_2,fmt="xg",label = 'Data with Taylor')

    plt.xlim(0,2.6)
    
    plt.legend()
    
    plt.xlabel('$\\beta$')
    
    plt.ylabel('$X$')
    
    ax.spines[['right', 'top']].set_visible(False)
    
    plt.savefig('ChiralResults/Plots/Susceptibility.svg')
    plt.show()
def plot_generic_data_beta(betas,N,SU, order, N_order, N_measure, N_thermal, observable_name, symbol):
    result = [0 for i in range(len(betas))]
    err = [0 for j in range(len(betas))]
    for i,beta in enumerate(betas):
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
        result[i],err[i] = Stats(np.load(file_name)).estimate()
    ax = plt.subplot(111)
    plt.xlabel('$\\beta$')
    
    plt.ylabel(symbol)
    ax.spines[['right', 'top']].set_visible(False)
    plt.xlim(0,max(betas)+0.1)
    ax.errorbar(x=betas,y=result,yerr= err,fmt="xg",label = 'Data with Taylor')
    plt.show()

def strong_3(beta):
    return 1/4 *(beta/2 + (7 *beta**3)/192 + (63 *beta**5)/20480)
def strong_2(beta):
    return 1/4 *((2 *beta)/3 + beta**2/6 + beta**3/27 + (25 *beta**4)/864 )
def strong_coupling(beta):
    return 1/2*beta+1/6*beta**3-1/6*beta**5
def strong_coupling_alt(beta):
    return 1/2*beta+1/6*beta**3+1/6*beta**5
def weak_3(beta):
    Q_1 = 0.0958876
    Q_2 = -0.067
    return 1 - 105/(512* beta**2) - 15/(16 *beta) - (705 + 3060* Q_1 + 36480 *Q_2)/(4096 *beta**3)
def weak_2(beta):
    Q_1 = 0.0958876
    Q_2 = -0.067
    return 1 - 7/(72* beta**2) - 2/(3* beta) - (137 + 396 *Q_2 + 684 *Q_1)/(2592 *beta**3)
def weak_coupling(beta):
    Q_1 = 0.0958876
    Q_2 = -0.067
    return 1- 3/(8*beta)*(1+1/(16*beta)+(1/64+3/16*Q_1+1/8*Q_2)*1/beta**2)