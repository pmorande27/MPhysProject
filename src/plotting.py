"""
Module used to make all hte relevants plots for the Chiral Model
"""
import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats
from scipy.optimize import curve_fit

def Greens_mom_2(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'ww corr'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values,values_err= np.load(file_name)
    xdata = np.arange(N+1)
    def model(t, dE):
        L = N
        return  (np.cosh(1/dE*(t-L/2)) -1) / (np.cosh(1/dE*L/2)-1)
    a =curve_fit(model,xdata,values/values[0],sigma=values_err/values[0])
    xs = [model(d,a[0][0]) for d in xdata]
    print(a[0][0])
    print(a)
    ax = plt.subplot(111)
    ax.plot(xdata,xs,label='Best Fit Model, $\epsilon=$'+str(round(a[0][0],3)))

    ax.errorbar(xdata,values/values[0],values_err/values[0],fmt='.k',label='HMC Data')
    plt.yscale('log')
    plt.xlabel('t')
    plt.ylabel('ww corr')
    plt.legend()
    ax.spines[['right', 'top']].set_visible(False)
    plt.title(" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal))
    file_name =  "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    plt.savefig(file_name)
    plt.show()

def plot_e_desinty(betas, model_params):
    """
    Function used to plot the Energy density for different
    parameters such as the lenght of the lattice or the group
    used
    """

    result = np.zeros(len(betas))

    error = np.zeros(len(betas))

    result_2 = np.zeros(len(betas))

    error_2 = np.zeros(len(betas))

    for i, beta in enumerate(betas):

        if model_params["su_parameter"] == 2:

            file_name = "ChiralResults/Processed/Action/Action beta = " + str(beta) \
                        + " N = " + str(model_params["n_length"])

            file_name += " SU = " + str(model_params["su_parameter"]) + \
                         " Order = " + str(0) + " N Order = " + str(0)

            file_name += " N measurements = " + str(model_params["n_measure"]) \
                        + " N Thermal = "

            file_name += str(model_params["n_thermal"])+'.npy'

            result[i],error[i] = np.load(file_name)


        file_name = "ChiralResults/Processed/Action/Action beta = " + \
                    str(betas[i]) + " N = " + str(model_params["n_length"]) + \
                     " SU = " + str(model_params["su_parameter"])+" Order = " \
                    + str(model_params["order"])+" N Order = " + str(model_params['n_order']) \
                    + " N measurements = " + str(model_params["n_measure"]) + " N Thermal = " \
                    + str(model_params["n_thermal"]) + '.npy'

        result_2[i], error_2[i] = np.load(file_name)

    axis = plt.subplot(111)

    if model_params["su_parameter"] == 2:

        axis.errorbar(x=betas, y=result, yerr=error,
                      fmt="xk", label='Data With Exact Exponetnial')

        strong_2s = [strong_coupling_alt(beta) for beta in np.linspace(0, 1, 100)]

        strong = [strong_coupling(beta) for beta in np.linspace(0, 1, 100)]

        plt.plot(np.linspace(0, 1, 100), strong_2s, 'b', label='Strong Coupling Alt')

        plt.plot(np.linspace(0, 1, 100), strong, 'b', label='Strong Coupling')

        weak = [weak_coupling(beta) for beta in np.linspace(0.5, 4.0, 1000)]

        plt.plot(np.linspace(0.5, 4.0, 1000), weak, 'r', label='Weak Coupling')

    axis.errorbar(x=betas, y=result_2, yerr=error_2, fmt="xg",
                  label='Data with Taylor of Order = ' + str(model_params["order"]) +
                  " and N_order = 10")

    if model_params["su_parameter"] == 3:

        values = [strong_2(beta) for beta in np.linspace(0, 2, 100)]

        plt.plot(np.linspace(0, 2, 100), values, 'b', label='Strong Coupling Expansion')

        values_w = [weak_2(beta) for beta in np.linspace(1, max(betas), 1000)]

        plt.plot(np.linspace(1, max(betas), 1000), values_w, 'r', label='Weak Coupling Expansion')

    if model_params["su_parameter"] == 4:

        values_w = [weak_3(beta) for beta in np.linspace(1, max(betas), 1000)]

        plt.plot(np.linspace(1, max(betas), 1000), values_w, 'r', label='Weak Coupling Expansion')

        strong = [strong_3(beta) for beta in np.linspace(0, 3, 100)]

        plt.plot(np.linspace(0, 3, 100), strong, 'b', label='Strong Coupling Expansion')

    plt.xlim(0, max(betas)+0.1)

    plt.legend()

    plt.xlabel('$\\beta$')

    plt.ylabel('$e$')

    axis.spines[['right', 'top']].set_visible(False)

    plt.savefig('ChiralResults/Processed/Plots/Energy_density_SU_' + str(model_params["su_parameter"]) \
                + "_N_" + str(model_params["n_length"]) + 'Order_' \
                + str(model_params["order"]) + 'N_order_' + str(model_params['n_order']) + '.svg')

    plt.show()

def plot_sus(betas, n_length, su_parameter):
    """
    Function used to plot the Susceptibility for different
    parameters such as the lenght of the lattice or the group
    used
    """


    result_2 = np.zeros(len(betas))

    error_2 = np.zeros(len(betas))

    for i, beta in enumerate(betas):

        file_name_2 = "ChiralResults/Processed/Susceptibility/Susceptibility 2" + \
            " beta = " + str(beta) +" N = "+str(n_length)+ " SU = " + str(su_parameter)+".npy"

        #values = np.load(file_name)

        values_2 = np.load(file_name_2)

        #result[i],error[i] = Stats(values).estimate()

        result_2[i], error_2[i] = values_2

    axis = plt.subplot(111)

    #ax.errorbar(x=betas,y=result,yerr= error,fmt="xk",label = 'Data')

    axis.errorbar(x=betas, y=result_2, yerr=error_2, fmt="xg", label='Data with Taylor')

    plt.xlim(0, 2.6)

    plt.legend()

    plt.xlabel('$\\beta$')

    plt.ylabel('$X$')

    axis.spines[['right', 'top']].set_visible(False)

    plt.savefig('ChiralResults/Processed/Plots/Susceptibility.svg')

    plt.show()

def plot_generic_data_beta(betas, model_params, observable_name, symbol):
    """
    Function used to plot generic data of an observable
    against beta for different parameters such as
    the lenght of the lattice or the group used
    """

    result = [0 for i in range(len(betas))]

    err = [0 for j in range(len(betas))]

    for i, beta in enumerate(betas):

        file_name = "ChiralResults/Processed" + observable_name + "/" + \
                    observable_name + " beta = " + str(beta) + " N = " + \
                    str(model_params["n_length"]) + \
                    " SU = " + str(model_params["su_parameter"]) + " Order = " \
                    + str(model_params["order"]) + " N Order = " \
                    + str(model_params["n_order"]) + " N measurements = " + \
                    str(model_params["n_measure"]) + \
                    " N Thermal = " + str(model_params["n_thermal"]) + '.npy'

        result[i], err[i] = np.load(file_name)

    axis = plt.subplot(111)

    plt.xlabel('$\\beta$')

    plt.ylabel(symbol)

    axis.spines[['right', 'top']].set_visible(False)

    plt.xlim(0, max(betas) + 0.1)

    axis.errorbar(x=betas, y=result, yerr=err, fmt="xg", label='Data with Taylor')

    plt.show()

def strong_3(beta):
    """
    Function for the strong couplig expansion of the
    internal energy density for the SU(4) group
    """

    return 1/4 * (beta/2 + (7 * beta**3)/192 + (63 * beta**5)/20480)

def strong_2(beta):
    """
    Function for the strong couplig expansion of the
    internal energy density for the SU(3) group
    """

    return 1/4 * ((2 * beta)/3 + beta**2/6 + beta**3/27 + (25 * beta**4)/864)

def strong_coupling(beta):
    """
    Function for the strong couplig expansion of the
    internal energy density for the SU(2) group
    """

    return 1/2 * beta+1/6 * beta**3-1/6 * beta**5

def strong_coupling_alt(beta):
    """
    Alternative Function for the strong couplig expansion
    of the internal energy density for the SU(3) group
    """

    return 1/2 * beta+1/6 * beta**3+1/6 * beta**5

def weak_3(beta):
    """
    Function for the weak couplig expansion of the
    internal energy density for the SU(4) group
    """

    q_one = 0.0958876

    q_two = -0.067

    return 1 - 105/(512 * beta**2) - 15/(16 * beta) - \
            (705 + 3060 * q_one + 36480 * q_two)/(4096 * beta**3)

def weak_2(beta):
    """
    Function for the weak couplig expansion of the
    internal energy density for the SU(3) group
    """

    q_one = 0.0958876

    q_two = -0.067

    return 1 - 7/(72 * beta**2) - 2/(3 * beta) - (137 + 396 * q_two + 684 * q_one)/(2592 * beta**3)

def weak_coupling(beta):
    """
    Function for the weak couplig expansion of the
    internal energy density for the SU(2) group
    """

    q_one = 0.0958876

    q_one = -0.067

    return 1- 3/(8 * beta) * (1 + 1/(16 * beta) + (1/64 + 3/16 * q_one + 1/8 * q_one) * 1/beta**2)

def plot_Greens_0_mom(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 0 Mom'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values,errors = np.load(file_name)
    xdata = np.arange(N+1)
    def model(t, dE):
        L = N
        return  (np.cosh(1/dE*(t-L/2)) -1) / (np.cosh(1/dE*L/2)-1)
    a =curve_fit(model,xdata,values/values[0])
    ys = [model(d,a[0][0]) for d in xdata]
    
    print(a)
    ax = plt.subplot(111)
    ax.errorbar(xdata,values,yerr=errors,fmt='.k')
    ax.plot(xdata,ys)
    ax.spines[['right', 'top']].set_visible(False)
    plt.show()
def plot_Greens_diags(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens Diags'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values,errors = np.load(file_name)
    xdata = np.arange(N+1)
    def model(t, dE):
        L = N
        return  (np.cosh(1/dE*(t-L/2)) -1) / (np.cosh(1/dE*L/2)-1)
    a =curve_fit(model,xdata,values/values[0])
    ys = [model(d,a[0][0]) for d in xdata]
    
    print(a)
    ax = plt.subplot(111)
    ax.errorbar(xdata,values,yerr=errors,fmt='.k')
    ax.plot(xdata,ys)
    ax.spines[['right', 'top']].set_visible(False)
    plt.show()
