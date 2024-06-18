"""
Module used to make all hte relevants plots for the Chiral Model
"""
import numpy as np
from matplotlib.ticker import MaxNLocator
import scipy as sci
import matplotlib.pyplot as plt
from Stats import Stats
from scipy.optimize import curve_fit

def correlation_length_automatic_fit(beta,N,SU,order,N_order,N_measure,N_thermal,lower,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], np.sqrt((ww_cor_err/ww_cor[0])**2 + (ww_cor/ww_cor[0]**2 * ww_cor_err[0])**2)
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)) / (np.cosh(N_2/xi) )
    cor_length,cor_length_err,reduced_chi2 = np.zeros(N_2),np.zeros(N_2),np.zeros(N_2)
    for i,upper_fit in enumerate(range(lower+1,N_2)):
    # perform the fit  
        one = ds <= upper_fit
        two = ds >=lower
        mask = one*two# fitting range
        popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True)
        cor_length[i] = popt[0] # in units of lattice spacing
        cor_length_err[i] = np.sqrt(pcov[0][0])

        r = cor[mask] - model(ds[mask], *popt)
        reduced_chi2[i] = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters
    i = np.argmin(abs(reduced_chi2-[1 for i in range(N_2)]))
    #print(range(1,N_2+1)[i])
    one = ds <= range(1,N_2+1)[i]
    two = ds >= lower
    mask = one*two# fitting range
    #print(mask)
    #mask = ds <= upper_fit
    reduced_chi2 = reduced_chi2[i]
    cor_length = cor_length[i]
    cor_length_err = cor_length_err[i]
    
    fig = plt.figure(figsize=(8,6))
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    axis.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black',label='Data')
    ds_fit = np.linspace(ds[mask][0], ds[mask][-1], 500)
    axis.plot(ds_fit, model(ds_fit,*popt),linestyle='dashed',linewidth=0.75, c='black', label='$\\xi = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    
    def model2(d,xi):
        return (np.cosh((d-N_2)/xi)-1) / (np.cosh(N_2/xi)-1 )
    cor_length,cor_length_err,reduced_chi2 = np.zeros(N_2),np.zeros(N_2),np.zeros(N_2)
    for i,upper_fit in enumerate(range(lower+1,N_2)):
    # perform the fit  
        one = ds <= upper_fit
        two = ds >=lower
        mask = one*two# fitting range
        popt, pcov = curve_fit(model2, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True)
        cor_length[i] = popt[0] # in units of lattice spacing
        cor_length_err[i] = np.sqrt(pcov[0][0])

        r = cor[mask] - model2(ds[mask], *popt)
        reduced_chi2[i] = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters
    i = np.argmin(abs(reduced_chi2-[1 for i in range(N_2)]))
    #print(range(1,N_2+1)[i])
    one = ds <= range(1,N_2+1)[i]
    two = ds >= lower
    mask = one*two# fitting range
    #print(mask)
    #mask = ds <= upper_fit
    reduced_chi2 = reduced_chi2[i]
    cor_length = cor_length[i]
    cor_length_err = cor_length_err[i]
    ds_fit = np.linspace(ds[mask][0], ds[mask][-1], 500)
    axis.plot(ds_fit, model2(ds_fit,*popt),linestyle='dashed',linewidth=0.75, c='red', label='$\\xi_{-1} = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    

    plt.yscale('log')
    plt.ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.svg'
    
    plt.savefig(file_name)
    plt.show()

def correlation_length_manual_fit(beta,N,SU,order,N_order,N_measure,N_thermal,lower,upper,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)) / (np.cosh(N_2/xi))
    upper_fit = upper
    lower_fit = lower
    one = ds <= upper_fit
    two = ds >=lower_fit
    mask = one*two# fitting range
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True, bounds=([0],[np.inf]))
    cor_length= popt[0] # in units of lattice spacing
    cor_length_err = np.sqrt(pcov[0][0])

    r = cor[mask] - model(ds[mask], *popt)
    reduced_chi2 = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters
    
    
    fig = plt.figure(figsize=(8,6))
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    axis.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    cor_length,cor_err = cor_length,cor_err
    ds_fit = np.linspace(ds[mask][0], ds[mask][-1], 500)
    axis.plot(ds_fit, model(ds_fit,*popt), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)-1) / (np.cosh(N_2/xi)-1)
    upper_fit = upper
    lower_fit = lower
    one = ds <= upper_fit
    two = ds >=lower_fit
    mask = one*two# fitting range
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True, bounds=([0],[np.inf]))
    cor_length= popt[0] # in units of lattice spacing
    cor_length_err = np.sqrt(pcov[0][0])

    r = cor[mask] - model(ds[mask], *popt)
    reduced_chi2 = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters
    axis.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    cor_length,cor_err = cor_length,cor_err
    ds_fit = np.linspace(ds[mask][0], ds[mask][-1], 500)
    axis.plot(ds_fit, model(ds_fit,*popt), c='red',linestyle='dashed',linewidth=0.75, label='$\\xi_{-1} = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    

    plt.yscale('log')
    plt.ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()

def plot_e_desinty(betas, model_params,accel = False):
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
            if accel == False:

                file_name = "ChiralResults/Processed/Action/Action beta = " + str(beta) \
                            + " N = " + str(model_params["n_length"])

                file_name += " SU = " + str(model_params["su_parameter"]) + \
                            " Order = " + str(0) + " N Order = " + str(0)

                file_name += " N measurements = " + str(model_params["n_measure"]) \
                            + " N Thermal = "

                file_name += str(model_params["n_thermal"])+'.npy'
            else:
                file_name = "ChiralResults/Processed/Action/Action beta = " + str(beta) \
                            + " N = " + str(model_params["n_length"])

                file_name += " SU = " + str(model_params["su_parameter"]) + \
                            " Order = " + str(0) + " N Order = " + str(0)

                file_name += " N measurements = " + str(model_params["n_measure"]) \
                            + " N Thermal = "

                file_name += str(model_params["n_thermal"])+" Accel.npy"
            result[i],error[i] = np.load(file_name)

        if accel == False:

            file_name = "ChiralResults/Processed/Action/Action beta = " + \
                    str(betas[i]) + " N = " + str(model_params["n_length"]) + \
                     " SU = " + str(model_params["su_parameter"])+" Order = " \
                    + str(model_params["order"])+" N Order = " + str(model_params['n_order']) \
                    + " N measurements = " + str(model_params["n_measure"]) + " N Thermal = " \
                    + str(model_params["n_thermal"]) + '.npy'
        else:
            file_name = "ChiralResults/Processed/Action/Action beta = " + \
                    str(betas[i]) + " N = " + str(model_params["n_length"]) + \
                     " SU = " + str(model_params["su_parameter"])+" Order = " \
                    + str(model_params["order"])+" N Order = " + str(model_params['n_order']) \
                    + " N measurements = " + str(model_params["n_measure"]) + " N Thermal = " \
                    + str(model_params["n_thermal"]) + " Accel.npy"

        result_2[i], error_2[i] = np.load(file_name)

    axis = plt.subplot(111)

    if model_params["su_parameter"] == 2:

        axis.errorbar(x=betas, y=result, yerr=error,
                      fmt="xk",color='black', label='Data With Exact Exponetnial')

        strong_2s = [strong_coupling_alt(beta) for beta in np.linspace(0, 1, 100)]

        strong = [strong_coupling(beta) for beta in np.linspace(0, 1, 100)]

        plt.plot(np.linspace(0, 1, 100), strong_2s, 'blue',linestyle='dashdot',linewidth=0.75, label='Strong Coupling Alt')

        plt.plot(np.linspace(0, 1, 100), strong, 'black',linestyle='--',linewidth=0.75, label='Strong Coupling')

        weak = [weak_coupling(beta) for beta in np.linspace(0.5, 4.0, 1000)]

        plt.plot(np.linspace(0.5, 4.0, 1000), weak, 'red',linestyle='dashed',linewidth=0.75, label='Weak Coupling')

    axis.errorbar(x=betas, y=result_2, yerr=error_2, fmt=".",color='black',
                  label='Data with Taylor of Order = ' + str(model_params["order"]) +
                  " and N_order = 10")

    if model_params["su_parameter"] == 3:

        values = [strong_2(beta) for beta in np.linspace(0, max(betas), 1000)]

        plt.plot(np.linspace(0, max(betas), 1000), values, 'b', label='Strong Coupling Expansion',linewidth=0.75,linestyle='dashed')

        values_w = [weak_2(beta) for beta in np.linspace(0.1, max(betas), 10000)]

        plt.plot(np.linspace(0.1, max(betas), 10000), values_w, 'r', label='Weak Coupling Expansion',linewidth=0.75,linestyle='dashed')

    if model_params["su_parameter"] == 4:

        values_w = [weak4(beta) for beta in np.linspace(0.1, max(betas), 10000)]

        plt.plot(np.linspace(0.1, max(betas), 10000), values_w, 'r', label='Weak Coupling Expansion',linewidth=0.75,linestyle='dashed')

        strong = [strong_3(beta) for beta in np.linspace(0, max(betas), 1000)]

        plt.plot(np.linspace(0, max(betas), 1000), strong, 'b', label='Strong Coupling Expansion',linewidth=0.75,linestyle='dashed')
    if model_params["su_parameter"] == 6:

        
        strong = [strong_6(beta) for beta in np.linspace(0, 8, 10000)]
        values_w = [weak6(beta) for beta in np.linspace(0.1, max(betas), 10000)]

        plt.plot(np.linspace(0.1, max(betas), 10000), values_w, 'r', label='Weak Coupling Expansion',linewidth=0.75,linestyle='dashed')

        plt.plot(np.linspace(0, 8, 10000), strong, 'b', label='Strong Coupling Expansion',linewidth=0.75,linestyle='dashed')
    if model_params["su_parameter"] == 9:

        
        strong = [strong_9(beta) for beta in np.linspace(0, 8, 10000)]
        values_w = [weak9(beta) for beta in np.linspace(0.1, max(betas), 10000)]

        plt.plot(np.linspace(0.1, max(betas), 10000), values_w, 'r', label='Weak Coupling Expansion',linewidth=0.75,linestyle='dashed')

        plt.plot(np.linspace(0, 8, 10000), strong, 'b', label='Strong Coupling Expansion',linewidth=0.75,linestyle='dashed')

    plt.xlim(0, max(betas) + 0.1)
    plt.ylim(0, 1.0)

    plt.legend()

    plt.xlabel('$\\beta$')

    plt.ylabel('$e$')

    axis.spines[['right', 'top']].set_visible(False)
    if accel == False:
        plt.savefig('ChiralResults/Processed/Plots/Energy_density_SU_' + str(model_params["su_parameter"]) \
                + "_N_" + str(model_params["n_length"]) + 'Order_' \
                + str(model_params["order"]) + 'N_order_' + str(model_params['n_order']) + '.svg')
    else:
        plt.savefig('ChiralResults/Processed/Plots/Energy_density_SU_' + str(model_params["su_parameter"]) \
                + "_N_" + str(model_params["n_length"]) + 'Order_' \
                + str(model_params["order"]) + 'N_order_' + str(model_params['n_order']) + ' Accel.svg')

    plt.show()



def plot_generic_data_beta(betas, model_params, observable_name, symbol, accel = False):
    """
    Function used to plot generic data of an observable
    against beta for different parameters such as
    the lenght of the lattice or the group used
    """

    result = [0 for i in range(len(betas))]

    err = [0 for j in range(len(betas))]

    for i, beta in enumerate(betas):
        if accel == False:

            file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                        observable_name + " beta = " + str(beta) + " N = " + \
                        str(model_params["n_length"]) + \
                        " SU = " + str(model_params["su_parameter"]) + " Order = " \
                        + str(model_params["order"]) + " N Order = " \
                        + str(model_params["n_order"]) + " N measurements = " + \
                        str(model_params["n_measure"]) + \
                        " N Thermal = " + str(model_params["n_thermal"]) + '.npy'
        else:
            file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                        observable_name + " beta = " + str(beta) + " N = " + \
                        str(model_params["n_length"]) + \
                        " SU = " + str(model_params["su_parameter"]) + " Order = " \
                        + str(model_params["order"]) + " N Order = " \
                        + str(model_params["n_order"]) + " N measurements = " + \
                        str(model_params["n_measure"]) + \
                        " N Thermal = " + str(model_params["n_thermal"]) + ' Accel.npy'

        result[i], err[i] = np.load(file_name)

    axis = plt.subplot(111)

    plt.xlabel('$\\beta$')

    plt.ylabel(symbol)

    axis.spines[['right', 'top']].set_visible(False)

    plt.xlim(0, max(betas) + 0.1)

    axis.errorbar(x=betas, y=result, yerr=err, fmt="xg", label='Data with Taylor')

    plt.show()
    plt.savefig("ChiralResults/Processed/Plots/Gen Obs"+".svg")
    

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
def strong_6(beta):
    return beta/12 +beta**3/864 +7*beta**5/103680
def strong_9(beta):
    b = beta/18
    return 1/4* (4 *b + 8 *b**3 + 24* b**5)
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

def weak6(beta):
    beta = beta/12
    q_one = 0.0958876
    q_two = -0.067
    a_one = 17/576
    n = 6
    a_two =  (3*n**4-14*n**2+20)/(768*n**4)+ q_one*(n**4-4*n**2+12)/(64*n**4) + q_two*(n**4-n**2+24)/(64*n**4)
    return 1-(n**2-1)/(8*n**2*beta)*(1 +a_one/beta + a_two/beta**2)  
def weak4(beta):
    beta = beta/8
    q_one = 0.0958876
    q_two = -0.067
    a_one = 7/256
    n = 4
    a_two =  (3*n**4-14*n**2+20)/(768*n**4)+ q_one*(n**4-4*n**2+12)/(64*n**4) + q_two*(n**4-n**2+24)/(64*n**4)
    return 1-(n**2-1)/(8*n**2*beta)*(1 +a_one/beta + a_two/beta**2)
def weak9(beta):
    beta = beta/18
    q_one = 0.0958876
    q_two = -0.067
    a_one = 7/256
    n = 9
    a_two =  (3*n**4-14*n**2+20)/(768*n**4)+ q_one*(n**4-4*n**2+12)/(64*n**4) + q_two*(n**4-n**2+24)/(64*n**4)
    return 1-(n**2-1)/(8*n**2*beta)*(1 +a_one/beta + a_two/beta**2)
def plot_Greens_0_mom(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens 0 mom'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    ww_cor,ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi) - 1) / (np.cosh(N_2/xi) - 1)

    # perform the fit  
    mask = cor > 0 # fitting range
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True)
    cor_length = popt[0] # in units of lattice spacing
    cor_length_err = np.sqrt(pcov[0][0])

    r = cor[mask] - model(ds[mask], *popt)
    reduced_chi2 = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters

    
    fig = plt.figure(figsize=(8,6))

    plt.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2)
    ds_fit = np.linspace(0, ds[mask][-1], 500)
    plt.plot(ds_fit, model(ds_fit,*popt), c='g', label='$\\xi = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    plt.yscale('log')
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    plt.show()
def plot_Greens_diags(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens diags'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    values,errors = np.load(file_name)
    ww_cor,ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi) - 1) / (np.cosh(N_2/xi) - 1)

    # perform the fit  
    mask = cor > -1 # fitting range
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True)
    cor_length = popt[0] # in units of lattice spacing
    cor_length_err = np.sqrt(pcov[0][0])

    r = cor[mask] - model(ds[mask], *popt)
    reduced_chi2 = np.sum((r/cor_err[mask])**2) / (mask.size - 1) # dof = number of observations - number of fitted parameters

    
    fig = plt.figure(figsize=(8,6))

    plt.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2)
    cor_length,cor_length_err = cor_length/np.sqrt(2),cor_length_err/np.sqrt(2)
    ds_fit = np.linspace(0, ds[mask][-1], 500)
    plt.plot(ds_fit, model(ds_fit,*popt), c='g', label='$\\xi = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    plt.yscale('log')
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    plt.show()
def plot_Greens_diags2(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens Diags'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    values,errors = np.load(file_name)
    xdata = np.arange(N-1)
    results =1/ np.abs(np.sqrt(2)*np.array([np.log(values[i+1])-np.log(values[i]) for i in range(N-1)]))
    print(results)
    
    ax = plt.subplot(111)
    #ax.errorbar(xdata,values,yerr=errors,fmt='.k')
    ax.plot(xdata,results)
    ax.spines[['right', 'top']].set_visible(False)
    plt.show()
def plot_iats(betas,N,SU,order,N_order,N_measure,N_thermal,corr_lens):
    observable_name = 'Susceptibility iat'
    l = len(corr_lens)
    vals_no_acc = np.zeros(l)
    vals_acc = np.zeros(l)
    err_no_acc = np.zeros(l)
    err_acc = np.zeros(l)
    for i,beta in enumerate(betas):
        file_name1 = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    
        file_name2 = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
        vals_no_acc[i],err_no_acc[i] = np.load(file_name1)
        vals_acc[i],err_acc[i] = np.load(file_name2)
    plt.errorbar(corr_lens,vals_no_acc,yerr=err_no_acc,fmt='^',capsize=2,label='No Accel',color='black')
    plt.errorbar(corr_lens,vals_acc,yerr=err_acc,fmt='.',capsize=2,label='Accel',color = 'black')
    log_a = np.log(corr_lens)
    log_iat = np.log(vals_acc)
    err_log = err_acc/vals_acc
    def linear(x, z, b):
        return z*x + b
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat,sigma=err_log, absolute_sigma=True)
    fitted_acc = [a**popt[0]*np.exp(popt[1]) for a in corr_lens]


    plt.plot(corr_lens,fitted_acc,label = 'Fit for accelration: $z = $' + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)),linestyle = 'dashed',color = 'blue',linewidth = 0.75)
    log_a = np.log(corr_lens)
    log_iat = np.log(vals_no_acc)
    err_log = err_no_acc/vals_no_acc
    def linear(x, z, b):
        return z*x + b
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat,sigma=err_log, absolute_sigma=True)
    fitted_acc = [a**popt[0]*np.exp(popt[1]) for a in corr_lens]
    plt.plot(corr_lens,fitted_acc,label = 'Fit for no accelration: $z = $' + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)),linestyle = 'dashed',color = 'black',linewidth = 0.75)

    vals2 = np.array([8.563446099663128308e+00, 8.948369438472782988e+00 ,9.919038263108260978e+00, 9.608982546817623316e+00 ,1.156136494640124646e+01, 1.259567673318285941e+01, 1.406568553927657028e+01])
    errs2 = np.array([3.259732171137674772e-01, 3.453898030471281855e-01, 4.032888832987974181e-01, 3.858292297862710440e-01, 5.090671684954585219e-01, 5.774897684578449431e-01, 6.814186249648421789e-01])
    log_iat = np.log([8.563446099663128308e+00, 8.948369438472782988e+00 ,9.919038263108260978e+00, 9.608982546817623316e+00 ,1.156136494640124646e+01, 1.259567673318285941e+01, 1.406568553927657028e+01])
    err_log = errs2/vals2
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat,sigma=err_log, absolute_sigma=True)
    fitted_acc = [a**popt[0]*np.exp(popt[1]) for a in corr_lens]
    plt.errorbar(corr_lens,vals2,yerr=errs2,fmt='x',capsize=2,label='Julian Data',color = 'black')
    plt.plot(corr_lens,fitted_acc,label = 'Fit for Julian accelration: $z = $' + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)),linestyle = 'dashed',color = 'red',linewidth = 0.75)



    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$\\xi [a]$')
    plt.ylabel('$\\tau_{iat}$')
    plt.savefig('ChiralResults/Processed/Plots/Susceptibility iat.svg')
    plt.show()
def comp_plot():
    N = 120
    l = 7.5
    m = 7
    A = 1
    B= 1
    def f(x):
        return A*(np.exp(-m*x)+np.exp(-m*(N-x))) + B*(np.exp(-l*x)+np.exp(-l*(N-x)))
    def m_e(x):
        return -np.log(f(x+1)/f(x))
    def m_c(x):
        return np.arccosh((f(x+1)+f(x-1))/(2*f(x)))
    xs = np.arange(N)
    plt.plot(xs,[m_e(x) for x in xs],label = 'Exp',color = 'black',linestyle = 'dashed',linewidth = 0.75)
    plt.plot(xs,[m_c(x) for x in xs],label = 'Cosh',color = 'red',linestyle = 'dashed',linewidth = 0.75)
    plt.legend()    
    plt.show()

def one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns):

    Probabilities = [np.exp(-chi_sq/2-k-n) for (chi_sq,k,n) in zip(chi_sqs,ks,ns)]
    Probabilities = np.array(Probabilities)
    Probabilities = Probabilities/np.sum(Probabilities)
    a = np.sum(Probabilities*a0s)
    sigmassq = np.sum(Probabilities*sigmassqs) +np.sum(Probabilities*a0s**2) - np.sum(Probabilities*a0s)**2 
    return a,sigmassq, Probabilities

def two_dimnesional_aic(a0s,a1s,sigmassqs,sigmasqss1,chi_sqs,ks,ns):

    Probabilities = [np.exp(-chi_sq/2-k-n) for (chi_sq,k,n) in zip(chi_sqs,ks,ns)]
    Probabilities = np.array(Probabilities)
    Probabilities = Probabilities/np.sum(Probabilities)
    a0 = np.sum(Probabilities*a0s)
    a1 = np.sum(Probabilities*a1s)

    sigmassq0 = np.sum(Probabilities*sigmassqs) +np.sum(Probabilities*a0s**2) - np.sum(Probabilities*a0s)**2 
    sigmassq1 = np.sum(Probabilities*sigmasqss1) +np.sum(Probabilities*a1s**2) - np.sum(Probabilities*a1s)**2
    return a0,sigmassq0,a1,sigmasqss1, Probabilities

def mass_plot_def(beta,N,SU,order,N_order,N_measure,N_thermal,lower_cs,upper_cs,lower_es,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    rel_err = cor_err / cor # relative error
    cor_1, cor_err_1 = np.roll(cor, -1), np.roll(cor_err, -1) # shift to d+1
    rel_err_1 = cor_err_1 / cor_1
    cor__1, cor_err__1 = np.roll(cor, 1), np.roll(cor_err, 1) # shift to d-1
    rel_err__1 = cor_err__1 / cor__1

    A, B = cor_1/cor, cor__1/cor
    x = (A+B)/2 
    m_eff = np.arccosh(x)
    print(len(m_eff))
    delta_x = 1/2 * (A*(rel_err_1 - rel_err) + B*(rel_err__1 - rel_err))
    # delta_x = A/2*(np.sqrt(rel_err_1**2 + rel_err**2)) + B/2*(np.sqrt(rel_err__1**2 + rel_err**2))
    m_eff_err = 1/np.sqrt(x**2-1) * delta_x

    cor_1 = np.roll(cor, -1) # shift to d+1
    m_eff = - np.log(cor_1 / cor)
    m_eff_err = np.roll(cor_err, -1)/cor_1 - cor_err/cor 
    #plt.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.', capsize=2)

    def model(x, m):
        return m
    ks = []
    a0s = []
    ns = []
    chi_sqs = []
    sigmassqs = []
    for j in upper_cs:
        for i in lower_cs:
            one = ds <= j
            two = ds >= i

            mask = one & two
            popt, pcov = curve_fit(model, ds[mask], m_eff[mask], sigma=m_eff_err[mask], absolute_sigma=True)
            a0s += [popt[0]]
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            ys = np.array([popt[0] for x in ds[mask]])
            r = m_eff[mask] - ys
            chi_sqs += [np.sum((r/m_eff_err[mask])**2)] 
    chi_sqs = np.array(chi_sqs)
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)

    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    print(m,np.sqrt(cov_m))
    

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    #axis.spines[['right', 'top']].set_visible(False)
    one = ds <= max(upper_cs)
    two = ds >= min(upper_cs)
    one2 = ds <= max(lower_cs)
    two2 = ds >= min(lower_cs)

    mask2 = one2 & two2

    mask = one & two
    axs[0].errorbar(ds, m_eff, yerr=m_eff_err, fmt='.',color='black', capsize=2,label='Cosh')
    axs[0].plot(ds, [m for s in ds], c='black',linestyle='dashed',linewidth=1.0,label='$\\xi^{cosh} = %.3f \pm %.3f$'%(1/m,  np.sqrt(cov_m)/m**2))
    axs[0].fill_between(ds, m - np.sqrt(cov_m), m + np.sqrt(cov_m), alpha=0.2, color='red')
    axs[0].fill_between(ds, 1, y2=-1,where= mask, alpha=0.2, color='red',label = 't maxs')
    axs[0].fill_between(ds, 1, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't mins cosh')
    lows = [lower_cs for i in range(len(upper_cs))]
    Probs = Probs.reshape(len(upper_cs),len(lower_cs))
    a0s = a0s.reshape(len(upper_cs),len(lower_cs))
    sigmassqs = sigmassqs.reshape(len(upper_cs),len(lower_cs))
    axs[0].legend()

    for i in range(len(upper_cs)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(upper_cs[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(upper_cs[i]))

    mask2 = one2 & two2
    axs[0].set_ylim(popt[0]-popt[0]/10, popt[0]+popt[0]/10)
    plt.ylim(popt[0]-popt[0]/100, popt[0]+popt[0]/100)
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.legend()
    plt.ylabel('m')
    plt.xlabel('t')
    plt.xlim(0,21+1)
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    plt.savefig(file_name)
    plt.show()

def correlation_length_def_fit(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    v, e = np.load(file_name)
    ww_cor, ww_cor_err = np.concatenate((v,[v[0]])),np.concatenate((e,[e[0]]))
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)) / (np.cosh(N_2/xi))
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=(24, 26))
            a0s += [popt[0]]
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = cor[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/cor_err[mask])**2)]
            print(np.sum((r/cor_err[mask])**2)/(mask.size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    xi, cov_xi, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    def err_model(d,xi,err_xi):
        return abs((1/np.cosh(N/xi)* ((-d + N)* np.sinh((d - N_2)/xi) + N_2* np.cosh((d - N_2)/xi)* np.tanh(N/xi)))/xi**2)*err_xi
    axs[0].plot(ds_fit, model(ds_fit,xi), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(xi, np.sqrt(cov_xi)))
    axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    axs[0].set_yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def Accel.svg'
    #axis.set_xticks(ds)

    plt.savefig(file_name)
    plt.show()

def correlation_length_def_fit_state_two(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)) / (np.cosh(N_2/xi))
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True)
            a0s += [popt[0]]
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = cor[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/cor_err[mask])**2)]
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    xi, cov_xi = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    axis.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    def err_model(d,xi,err_xi):
        return abs((1/np.cosh(N/xi)* ((-d + N)* np.sinh((d - N_2)/xi) + N_2* np.cosh((d - N_2)/xi)* np.tanh(N/xi)))/xi**2)*err_xi
    axis.plot(ds_fit, model(ds_fit,xi), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(xi, np.sqrt(cov_xi)))
    axis.fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axis.fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axis.fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    

    plt.ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def Accel.svg'
    #axis.set_xticks(ds)

    plt.savefig(file_name)
    plt.show()

def plot_exponential(beta, N, SU, order, N_order, N_measure,N_thermal, Hot_start = True,accel =False):
    observable_name = "Exponential H"

    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    N_taus,averages,errors,rates,averages_two,errors_two = np.load(file_name)  
    delta_t = 1/np.array(N_taus)
    fig = plt.figure()
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    def model(x,A):
        return A*x**4
    popt, pcov = curve_fit(model, delta_t, averages_two, sigma=errors, absolute_sigma=True)
    r = averages_two - model(delta_t, *popt)
    chi_sqs = np.sum((r/errors_two)**2)/(len(delta_t)-1)


    axs[0].errorbar(delta_t,averages,yerr=errors,fmt='.',capsize=2,color = 'black')
    axs[0].plot([0,100],[1,1],linestyle='dashed',color='black')
    axs[0].set_ylabel('$\\langle \exp\{-H\}  \\rangle$')
    axs[1].errorbar(delta_t,averages_two,yerr=errors_two,fmt='.',capsize=2,color = 'black')
    t = np.linspace(0,max(delta_t),1000)
    axs[1].plot(t,model(t,*popt),linestyle='dashed',color='red',label='$x^4$ fit \n $\\frac{\chi_{sq}}{DOF} =%.3f $'%(chi_sqs))
    axs[1].legend()
    axs[1].set_ylabel('$\\Delta H$')
    axs[2].plot(delta_t, rates, 'o',color = 'black')
    axs[2].set_ylabel('$P_{acc}$')
    dd = min(delta_t)/10
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[1].spines[['right']].set_visible(False)

    axs[2].spines[['right']].set_visible(False)
    dd2 = max(delta_t)/10
    plt.xlim(min(delta_t)-dd,max(delta_t)+dd2)
    plt.xlabel('$\\Delta \\tau$')
    if accel == False:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+ " beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    
    plt.savefig(file_name)
    plt.show()

def plot_accH(beta, N, SU, order, N_order, N_measure,N_thermal, Hot_start = True,accel =False):
    observable_name = "Exponential H"

    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    N_taus,averages,errors,rates,averages_two,errors_two = np.load(file_name)  
    plt.errorbar(averages_two,rates,xerr=errors_two,fmt='.',capsize=2,color = 'black')
    ts = np.linspace(0,max(averages_two),1000,endpoint=True)  
    plt.plot(ts,sci.special.erfc(np.sqrt(ts)/2),linestyle='dashed',color='black',label = '$erfc(\\frac{\sqrt {\\langle \\Delta H \\rangle}}{2})$')
    plt.xlabel('$\\langle \\Delta H \\rangle$')
    plt.ylabel('$P_{acc}$')
    plt.legend()
    observable_name = "Acceptance Probability Delta H"
    if accel == False:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+ " beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    
    plt.savefig(file_name)
    plt.show()
from Stats import Stats


def mass_plot_def_2(beta,N,SU,order,N_order,N_measure,N_thermal,lower_cs,upper_cs,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    N_2 = int(N/2) 
    ds = np.arange(N_2+1)
    vals= np.load(file_name)
    Stat = Stats(vals)
    def function_3(x):
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(x, -1)
        x_2 = np.roll(x, 1)
        return - np.log(x_1/x)
    def function_2(x):
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)
        return - np.log(x_1/cor)
    def function(x):
        #X = np.concatenate((x,[x[0]]))
        #cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(x, -1)
        x_2 = np.roll(x, 1)
        return np.arccosh((x_1+x_2) / (2*x))
    def function_4(x):
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)


        return np.arccosh((x_1+x_2) / (2*cor))
    m_eff, m_eff_err,cov= Stat.jacknife_arbitrary_function_2D_2_2(N,100, function_4)
    print(cov)
    #plt.errorbar(np.arange(N), m_eff, yerr=m_eff_err, fmt='.',color='red', capsize=2,label='Cosh')
    #m_eff = np.concatenate((m_eff,[m_eff[0]]))
    #m_eff_err = np.concatenate((m_eff_err,[m_eff_err[0]]))
    #m_eff = 1/2 * (m_eff[:N_2+1] +abs( m_eff[N_2:][::-1]))
    #m_eff_err = np.sqrt(m_eff_err[:N_2+1]**2 + m_eff_err[N_2:][::-1]**2) / np.sqrt(2)
    #print(np.arange(N))
    #print(m_eff)
    plt.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.',color='black', capsize=2,label='1')
    
    m_eff, m_eff_err= Stat.jacknife_arbitrary_function_2D_2(N,100, function)
    
    #plt.errorbar(np.arange(N), m_eff, yerr=m_eff_err, fmt='.',color='red', capsize=2,label='Cosh')
    m_eff = np.concatenate((m_eff,[m_eff[0]]))
    m_eff_err = np.concatenate((m_eff_err,[m_eff_err[0]]))
    m_eff = 1/2 * (m_eff[:N_2+1] +abs( m_eff[N_2:][::-1]))
    m_eff_err = np.sqrt(m_eff_err[:N_2+1]**2 + m_eff_err[N_2:][::-1]**2) / np.sqrt(2)
    plt.errorbar(ds+0.01, m_eff, yerr=m_eff_err, fmt='.',color='red', capsize=2,label='2')

    
    def model(x, m):
        return m
    ks = []
    a0s = []
    ns = []
    chi_sqs = []
    sigmassqs = []
    
    one = ds <= upper_cs
    two = ds >= lower_cs

    mask = one & two
           
    #popt, pcov = curve_fit(model, ds[mask], m_eff[mask], sigma=m_eff_err[mask], absolute_sigma=True)
            
   

    #m, cov_m = popt[0],pcov[0][0]
    #print(m,np.sqrt(cov_m))
    

   
    #axis.spines[['right', 'top']].set_visible(False)
   

    #plt.plot(ds[mask], [m for s in ds[mask]], c='black',linestyle='dashed',linewidth=1.0,label='$\\xi^{cosh} = %.3f \pm %.3f$'%(1/m,  np.sqrt(cov_m)/m**2))
    #print(np.sqrt(cov_m)/m**2)
    plt.legend()
    #plt.ylim(popt[0]-popt[0]/10, popt[0]+popt[0]/10)
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.legend()
    plt.ylabel('m')
    plt.xlabel('t')
    plt.ylim(0,0.1)
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    #ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    ww_cor_1 = np.roll(ww_cor, -1)
    ww_cor_err_1 = np.roll(ww_cor_err, -1)
    m_eff = -np.log(ww_cor_1/ww_cor)
    m_eff_err = ww_cor_err_1/ww_cor_1 - (ww_cor_err/ww_cor)
    xdata = np.arange(N+1)
    m_eff[len(m_eff)-1] = m_eff[0]
    m_eff_err[len(m_eff_err)-1] = m_eff_err[0]

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) +0.2# wall separations covering half the lattice length
    m_eff = 1/2 * (m_eff[:N_2+1] +abs( m_eff[N_2:][::-1]))
    m_eff_err = np.sqrt(m_eff_err[:N_2+1]**2 + m_eff_err[N_2:][::-1]**2) / np.sqrt(2)

    
    #plt.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.', capsize=2)
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    plt.xlim(0,20)
    plt.show()

def cor_plot_def_2(beta,N,SU,order,N_order,N_measure,N_thermal,lower_cs,upper_cs,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    N_2 = int(N/2) 
    ds = np.arange(N_2+1)
    vals= np.load(file_name)
    Stat = Stats(vals)
    print('a')
    def function(x):
        
       
        

        return x
    ww_cor,ww_cor_err = np.zeros(N),np.zeros(N)
    for j in range(N):
        Stat = Stats(vals[j])
        ww_cor[j],ww_cor_err[j] = Stat.jacknife_arbitrary(500)
        print(j)
    ww_cor, ww_cor_err = np.concatenate((ww_cor,[ww_cor[0]])), np.concatenate((ww_cor_err,[ww_cor_err[0]]))
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    plt.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='blue',label='jacknife 1')
    #Stat = Stats(vals)
    #ww_cor, ww_cor_err= Stat.jacknife_function_2D(N, function)
    #ww_cor, ww_cor_err = np.concatenate((ww_cor,[ww_cor[0]])), np.concatenate((ww_cor_err,[ww_cor_err[0]]))
    #ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    #xdata = np.arange(N+1)

    """N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    plt.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black',label='jacknife')"""
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    

    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    plt.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='red',label='no jacknife')
    plt.legend()
    plt.xlim(0,20)
    plt.yscale('log')
    plt.show()

def mass_plot_def_jack(beta,N,SU,order,N_order,N_measure,N_thermal,lower_cs,upper_cs,lower_es,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    N_2 = int(N/2) 
    ds = np.arange(N_2+1)
    vals= np.load(file_name)
    Stat = Stats(vals)
    def function_2(x):
        #X = np.concatenate((x,[x[0]]))
        #cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(x, -1)
        x_2 = np.roll(x, 1)
        return - np.log(x_1/x)
    def function(x):
        #X = np.concatenate((x,[x[0]]))
        #cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(x, -1)
        x_2 = np.roll(x, 1)

        return np.arccosh((x_1+x_2) / (2*x))
    m_eff, m_eff_err= Stat.jacknife_arbitrary_function_2D(N,100, function)
    #plt.errorbar(np.arange(N), m_eff, yerr=m_eff_err, fmt='.',color='red', capsize=2,label='Cosh')
    m_eff = np.concatenate((m_eff,[m_eff[0]]))
    m_eff_err = np.concatenate((m_eff_err,[m_eff_err[0]]))
    m_eff = 1/2 * (m_eff[:N_2+1] +abs( m_eff[N_2:][::-1]))
    m_eff_err = np.sqrt(m_eff_err[:N_2+1]**2 + m_eff_err[N_2:][::-1]**2) / np.sqrt(2)
    
    #plt.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.',color='black', capsize=2,label='Cosh')
    def model(x, m):
        return m
    ks = []
    a0s = []
    ns = []
    chi_sqs = []
    sigmassqs = []
    for j in upper_cs:
        for i in lower_cs:
            one = ds <= j
            two = ds >= i

            mask = one & two
            popt, pcov = curve_fit(model, ds[mask], m_eff[mask], sigma=m_eff_err[mask], absolute_sigma=True)
            a0s += [popt[0]]
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            ys = np.array([popt[0] for x in ds[mask]])
            r = m_eff[mask] - ys
            chi_sqs += [np.sum((r/m_eff_err[mask])**2)] 
    chi_sqs = np.array(chi_sqs)
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)

    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    print(m,np.sqrt(cov_m))
    

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    #axis.spines[['right', 'top']].set_visible(False)
    one = ds <= max(upper_cs)
    two = ds >= min(upper_cs)
    one2 = ds <= max(lower_cs)
    two2 = ds >= min(lower_cs)

    mask2 = one2 & two2

    mask = one & two
    axs[0].errorbar(ds, m_eff, yerr=m_eff_err, fmt='.',color='black', capsize=2,label='Cosh')
    axs[0].plot(ds, [m for s in ds], c='black',linestyle='dashed',linewidth=1.0,label='$\\xi^{cosh} = %.3f \pm %.3f$'%(1/m,  np.sqrt(cov_m)/m**2))
    axs[0].fill_between(ds, m - np.sqrt(cov_m), m + np.sqrt(cov_m), alpha=0.2, color='red')
    axs[0].fill_between(ds, 1, y2=-1,where= mask, alpha=0.2, color='red',label = 't maxs')
    axs[0].fill_between(ds, 1, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't mins cosh')
    lows = [lower_cs for i in range(len(upper_cs))]
    Probs = Probs.reshape(len(upper_cs),len(lower_cs))
    a0s = a0s.reshape(len(upper_cs),len(lower_cs))
    sigmassqs = sigmassqs.reshape(len(upper_cs),len(lower_cs))
    axs[0].legend()

    for i in range(len(upper_cs)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(upper_cs[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(upper_cs[i]))

    mask2 = one2 & two2
    axs[0].set_ylim(popt[0]-popt[0]/10, popt[0]+popt[0]/10)
    plt.ylim(popt[0]-popt[0]/100, popt[0]+popt[0]/100)
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.legend()
    plt.ylabel('m')
    plt.xlabel('t')
    plt.xlim(0,21+1)
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    #plt.savefig(file_name)
    plt.show()



def correlation_length_def_fit_Julian(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    ww_cor, ww_cor_err = np.load(file_name)
    ww_cor, ww_cor_err = ww_cor/ww_cor[0], ww_cor_err/ww_cor[0]
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    # store processed correlation function data
    def model(d,xi):
        return (np.cosh((d-N_2)/xi)) / (np.cosh(N_2/xi))
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=(0, 100))
            a0s += [popt[0]]
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = cor[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/cor_err[mask])**2)]
            print([np.sum((r/cor_err[mask])**2)][0]/(cor[mask].size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    xi, cov_xi, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    def err_model(d,xi,err_xi):
        return abs((1/np.cosh(N/xi)* ((-d + N)* np.sinh((d - N_2)/xi) + N_2* np.cosh((d - N_2)/xi)* np.tanh(N/xi)))/xi**2)*err_xi
    axs[0].plot(ds_fit, model(ds_fit,xi), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(xi, np.sqrt(cov_xi)))
    axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    axs[0].set_yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' def Accel.svg'
    #axis.set_xticks(ds)

    plt.savefig(file_name)
    plt.show()


def plot_log_mass_jack(beta,N,SU,order,N_order,N_measure,N_thermal,t_max,t_min,accel = False):
    observable_name = 'log mass'
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
    mes  = np.load(file_name)
    mass, mass_err = mes[0], mes[1]
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'cov.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel cov.npy'
    
    cov = np.load(file_name)
    plt.errorbar(np.arange(int(N/2)+1), mass, yerr=mass_err, fmt='.',color='black', capsize=2)
    ds = np.arange(int(N/2)+1)
    def model(x,m):
        return m
    
    one = ds < t_max
    two = ds >= t_min
    mask = one & two
    cov_mask = np.zeros((t_max-t_min,t_max-t_min))

    for i in range(t_max-t_min):
        for j in range(t_max-t_min):
            print(i+t_min,j+t_min)
            cov_mask[i,j] = cov[i+t_min,j+t_min]
    m, pcov = curve_fit(model, np.arange(int(N/2)+1)[mask], mass[mask], sigma=mass_err[mask], absolute_sigma=True)
    plt.plot(ds, [m[0] for s in ds], c='black',linestyle='dashed',linewidth=1.0)
    
    print(1/m[0],np.sqrt(pcov[0,0])/m[0]**2)
    plt.ylim(m[0]-m[0]/10,m[0]+m[0]/10)
    plt.show()
def plot_cosh_mass_jack(beta,N,SU,order,N_order,N_measure,N_thermal,t_max,t_min,accel = False):
    observable_name = 'cosh mass'
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
    mass, mass_err = np.load(file_name)
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'cov.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel cov.npy'
    cov = np.load(file_name)
    plt.errorbar(np.arange(int(N/2)+1), mass, yerr=mass_err, fmt='.',color='black', capsize=2)
    ds = np.arange(int(N/2)+1)
    def model(x,m):
        return m
    
    one = ds < t_max
    two = ds >= t_min

    mask = one & two
    cov_mask = np.zeros((t_max-t_min,t_max-t_min))
    
    for i in range(t_max-t_min):
        for j in range(t_max-t_min):
            cov_mask[i,j] = cov[i+t_min,j+t_min]
    m, pcov = curve_fit(model, np.arange(int(N/2)+1)[mask], mass[mask], sigma=mass_err[mask], absolute_sigma=True)
    chisq = 0
    for i in range(t_max-t_min):
        chisq += (mass[i+t_min] - m[0])*(mass[i+t_min] - m[0])/mass_err[i+t_min]**2
    print(m,np.sqrt(pcov[0,0]))
    plt.plot(ds, [m[0] for s in ds], c='black',linestyle='dashed',linewidth=1.0)
    
    plt.ylim(m[0]-m[0]/10,m[0]+m[0]/10)
    plt.show()
def correlation_length_fit_two_params(beta,N,SU,order,N_order,N_measure,N_thermal,lower,upper,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    v,e = np.load(file_name)
    ww_cor, ww_cor_err = np.concatenate((v,[v[0]])),np.concatenate((e,[e[0]]))

    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi,A):
        return A*(np.exp(-d/xi) + np.exp(-(N-d)/xi))
    upper_fit = upper
    lower_fit = lower
    one = ds <= upper_fit
    two = ds >=lower_fit
    mask = one*two# fitting range
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=[(0, 0), (100, np.Infinity)])
    cor_length= popt[0] # in units of lattice spacing
    A = popt[1]
    cor_length_err = np.sqrt(pcov[0][0])

    r = (cor[mask] - model(ds[mask], *popt))/cor_err[mask]
    reduced_chi2 = np.sum(r**2) / (ds[mask].size - 2) # dof = number of observations - number of fitted parameters
    
    fig = plt.figure(figsize=(8,6))
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    axis.errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    cor_length,cor_err = cor_length,cor_err
    ds_fit = np.linspace(ds[mask][0], ds[mask][-1], 500)
    print(A)

    axis.plot(ds_fit, model(ds_fit,*popt), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_length, cor_length_err, reduced_chi2))
    
    

    plt.yscale('log')
    plt.ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.ylabel('wall wall correlation $C_{ww}(d)$')
    plt.legend(prop={'size':12}, frameon=True, loc='upper right')
    observable_name = 'correlation length'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(cor_length,cor_length_err))
    observable_name = 'ww corr'
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()

def correlation_length_fit_two_params_2(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    v, e = np.load(file_name)
    ww_cor, ww_cor_err = np.concatenate((v,[v[0]])),np.concatenate((e,[e[0]]))
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi,A):
        return A*(np.exp(-d/xi) + np.exp(-(N-d)/xi))
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    a1s = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=(0, 100))
            a0s += [popt[0]]
            a1s = [popt[1]]
            sigmassqs += [pcov[0][0]]
            sigmassqs1 += [pcov[1][1]]
            ks += [2]
            ns += [len(ds) - len(ds[mask])]
            r = cor[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/cor_err[mask])**2)]
            print(np.sum((r/cor_err[mask])**2)/(cor[mask].size-2))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    a1s = np.array(a1s)
    sigmassqs1 = np.array(sigmassqs1)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    xi, cov_xi,A,cov_A, Probs = two_dimnesional_aic(a0s,a1s,sigmassqs,sigmassqs1,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    def err_model(d,xi,err_xi):
        return abs((1/np.cosh(N/xi)* ((-d + N)* np.sinh((d - N_2)/xi) + N_2* np.cosh((d - N_2)/xi)* np.tanh(N/xi)))/xi**2)*err_xi
    axs[0].plot(ds_fit, model(ds_fit,xi,A), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(xi, np.sqrt(cov_xi)))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 100, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 100, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(10**(-4),10**2)
    plt.xlabel(r'wall separation $d$ [$a$]')
    axs[0].set_yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

   

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()


def correlation_length_fit_two_params_2_state_2(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    v, e = np.load(file_name)
    ww_cor, ww_cor_err = np.concatenate((v,[v[0]])),np.concatenate((e,[e[0]]))
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,xi,A):
        return A*(np.exp(-d/xi) + np.exp(-(N-d)/xi))
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    a1s = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=(0, np.infty ))
            a0s += [popt[0]]
            a1s = [popt[1]]
            
            sigmassqs += [pcov[0][0]]
            sigmassqs1 += [pcov[1][1]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = cor[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/cor_err[mask])**2)]
            print(np.sum((r/cor_err[mask])**2)/(ds[mask].size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    a1s = np.array(a1s)
    sigmassqs1 = np.array(sigmassqs1)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    xi, cov_xi,A,cov_A, Probs = two_dimnesional_aic(a0s,a1s,sigmassqs,sigmassqs1,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    def err_model(d,xi,err_xi):
        return abs((1/np.cosh(N/xi)* ((-d + N)* np.sinh((d - N_2)/xi) + N_2* np.cosh((d - N_2)/xi)* np.tanh(N/xi)))/xi**2)*err_xi
    axs[0].plot(ds_fit, model(ds_fit,xi,A), c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(xi, np.sqrt(cov_xi)))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(10**(-4),10**0)
    plt.xlabel(r'wall separation $d$ [$a$]')
    axs[0].set_yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

   

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    m, m_err = 1/xi,np.sqrt(cov_xi)/xi**2
    print(1/xi,np.sqrt(cov_xi)/xi**2)
    

    plt.savefig(file_name)
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    print(m)
    np.save(file_name,(m ,m_err))
    plt.show()

def correlation_length_fit_two_params_2_state_2_2(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    v, e = np.load(file_name)
    ww_cor, ww_cor_err = np.concatenate((v,[v[0]])),np.concatenate((e,[e[0]]))
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    cor = 1/2 * (ww_cor[:N_2+1] + ww_cor[N_2:][::-1])
    cor_err = np.sqrt(ww_cor_err[:N_2+1]**2 + ww_cor_err[N_2::-1]**2) / np.sqrt(2)

    # store processed correlation function data
    def model(d,m,A):
        return A*(np.exp(-d*m) + np.exp(-(N-d)*m))
   
    
    
    one = ds <= uppers
    two = ds >= lowers
    mask = one*two
    popt, pcov = curve_fit(model, ds[mask], cor[mask], sigma=cor_err[mask], absolute_sigma=True,bounds=[(0,0),(np.inf,np.inf)] )
    a0s = popt[0]
    a1s = popt[1]
            
    sigmassqs = pcov[0][0]
    sigmassqs1 = pcov[1][1]
           
    r = cor[mask] - model(ds[mask], *popt)
    chi_sqs = np.sum((r/cor_err[mask])**2)
    print(np.sum((r/cor_err[mask])**2)/(cor_err[mask].size-2))
    #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
   
    m, cov_m,A,cov_A,  = a0s,sigmassqs,a1s,sigmassqs1

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, cor, yerr=cor_err, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 1000,endpoint=True)
    one = ds_fit <= uppers
    two = ds_fit >= lowers
    mask = one*two
    


    mask = one & two
    axs[0].plot(ds_fit[mask], model(ds_fit[mask],m,A), c='black',linestyle='dashed',linewidth=0.75, label='$m = %.3f \pm %.3f$'%(m, np.sqrt(cov_m)))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    
   
    

    axs[0].set_ylim(10**(-4),10**3)
    plt.xlabel(r'wall separation $d$ [$a$]')
    axs[0].set_yscale('log')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')

   

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()

def plot_log_mass_jack(beta,N,SU,order,N_order,N_measure,N_thermal,t_max,t_min,accel = False):
    observable_name = 'log mass'
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
    mes  = np.load(file_name)
    mass, mass_err = mes[0], mes[1]
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'cov.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel cov.npy'
    
    cov = np.load(file_name)
    plt.errorbar(np.arange(int(N/2)+1), mass, yerr=mass_err, fmt='.',color='black', capsize=2)
    ds = np.arange(int(N/2)+1)
    def model(x,m):
        return m
    
    one = ds < t_max
    two = ds >= t_min
    mask = one & two
    cov_mask = np.zeros((t_max-t_min,t_max-t_min))

    for i in range(t_max-t_min):
        for j in range(t_max-t_min):
            print(i+t_min,j+t_min)
            cov_mask[i,j] = cov[i+t_min,j+t_min]
    m, pcov = curve_fit(model, np.arange(int(N/2)+1)[mask], mass[mask], sigma=mass_err[mask], absolute_sigma=True)
    plt.plot(ds, [m[0] for s in ds], c='black',linestyle='dashed',linewidth=1.0)
    
    print(1/m[0],np.sqrt(pcov[0,0])/m[0]**2)
    plt.ylim(m[0]-m[0]/10,m[0]+m[0]/10)
    plt.show()
def plot_cosh_mass_jack_state_2(beta,N,SU,order,N_order,N_measure,N_thermal,t_max,t_min,accel = False):
    observable_name = 'cosh mass state 2'
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
    mass, mass_err = np.load(file_name)
    if accel == False:

        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'cov.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+'/'+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel cov.npy'
    cov = np.load(file_name)
    plt.errorbar(np.arange(int(N/2)+1), mass, yerr=mass_err, fmt='.',color='black', capsize=2)
    ds = np.arange(int(N/2)+1)
    def model(x,m):
        return m
    
    one = ds < t_max
    two = ds >= t_min

    mask = one & two
    cov_mask = np.zeros((t_max-t_min,t_max-t_min))
    
    for i in range(t_max-t_min):
        for j in range(t_max-t_min):
            cov_mask[i,j] = cov[i+t_min,j+t_min]
    m, pcov = curve_fit(model, np.arange(int(N/2)+1)[mask], mass[mask], sigma=mass_err[mask], absolute_sigma=True)
    chisq = 0
    for i in range(t_max-t_min):
        chisq += (mass[i+t_min] - m[0])*(mass[i+t_min] - m[0])/mass_err[i+t_min]**2
    print(chisq/(t_max-t_min-1))
    print(m,np.sqrt(pcov[0,0]))

    plt.plot(ds, [m[0] for s in ds], c='black',linestyle='dashed',linewidth=1.0)
    
    plt.ylim(m[0]-m[0]/10,m[0]+m[0]/10)
    plt.show()

def plot_cosh_mass_jack_def(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'cosh mass'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    mass, e = np.load(file_name)
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length

    # store processed correlation function data
    def model(d,m):
        return m
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], mass[mask], sigma=e[mask], absolute_sigma=True,bounds=(0, np.infty ))
            a0s += [popt[0]]
            
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = mass[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/e[mask])**2)]
            print(np.sum((r/e[mask])**2)/(mass[mask].size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, mass, yerr=e, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    axs[0].plot(ds_fit,[m for i in ds_fit], c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(1/m, np.sqrt(cov_m)/m**2))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(m-m,m+m)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

    result = (1/m, np.sqrt(cov_m)/m**2)
    observable_name = 'correlation length'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,result)
    print(m)
    result = (m, np.sqrt(cov_m))
    observable_name = 'mass state 1'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,result)
    observable_name = 'cosh mass'

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()
def plot_log_mass_jack_def(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'log mass'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    mass, e = np.load(file_name)
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length

    # store processed correlation function data
    def model(d,m):
        return m
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], mass[mask], sigma=e[mask], absolute_sigma=True,bounds=(0, 3 ))
            a0s += [popt[0]]
            
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = mass[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/e[mask])**2)]
            print(np.sum((r/e[mask])**2)/(mass[mask].size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, mass, yerr=e, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    axs[0].plot(ds_fit,[m for i in ds_fit], c='black',linestyle='dashed',linewidth=0.75, label='$\\xi = %.3f \pm %.3f$'%(1/m, np.sqrt(cov_m)/m**2))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(m-m,m+m)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

   

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()

def plot_cosh_mass_jack_def_state_2(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'cosh mass state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    mass, e = np.load(file_name)
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length
    print(mass,e)
    # store processed correlation function data
    def model(d,m):
        return m
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], mass[mask], sigma=e[mask], absolute_sigma=True,bounds=(0, 3 ))
            a0s += [popt[0]]
            
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = mass[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/e[mask])**2)]
            print(np.sum((r/e[mask])**2)/(mass[mask].size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, mass, yerr=e, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    axs[0].plot(ds_fit,[m for i in ds_fit], c='black',linestyle='dashed',linewidth=0.75, label='$m = %.3f \pm %.3f$'%(m, np.sqrt(cov_m)))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(m-m/10,m+m/10)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')
    observable_name = 'mass state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    print(m)
    np.save(file_name,(m, np.sqrt(cov_m)))
    observable_name = 'cosh mass state 2'

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()
def plotlog_log_mass_jack_def_state_2(beta,N,SU,order,N_order,N_measure,N_thermal,lowers,uppers,accel = False):
    observable_name = 'log mass state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    mass, e = np.load(file_name)
    
    xdata = np.arange(N+1)

    N_2 = int(N/2) 
    ds = np.arange(N_2+1) # wall separations covering half the lattice length

    # store processed correlation function data
    def model(d,m):
        return m
    a0s = []
    sigmassqs = []
    ks = []
    ns = []
    chi_sqs = []
    sigmassqs1 = []
    
    for i in uppers:
        for j in lowers:
            one = ds <= i
            two = ds >= j
            mask = one*two
            popt, pcov = curve_fit(model, ds[mask], mass[mask], sigma=e[mask], absolute_sigma=True,bounds=(0, 3 ))
            a0s += [popt[0]]
            
            sigmassqs += [pcov[0][0]]
            ks += [1]
            ns += [len(ds) - len(ds[mask])]
            r = mass[mask] - model(ds[mask], *popt)
            chi_sqs += [np.sum((r/e[mask])**2)]
            print(np.sum((r/e[mask])**2)/(mass[mask].size-1))
            #print([np.sum((r/cor_err[mask])**2)][0]/(mask.size-1))
    a0s = np.array(a0s)
    sigmassqs = np.array(sigmassqs)
    ks = np.array(ks)
    ns = np.array(ns)
    chi_sqs = np.array(chi_sqs)
    m, cov_m, Probs = one_dimnesional_aic(a0s,sigmassqs,chi_sqs,ks,ns)
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    axs[0].spines[['right', 'top']].set_visible(False)
    axs[0].errorbar(ds, mass, yerr=e, fmt='.', capsize=2,color='black')
    ds_fit = np.linspace(ds[0], ds[-1], 500)
    
    one = ds <= max(uppers)
    two = ds >= min(uppers)
    one2 = ds <= max(lowers)
    two2 = ds >= min(lowers)

    mask2 = one2 & two2

    mask = one & two
    axs[0].plot(ds_fit,[m for i in ds_fit], c='black',linestyle='dashed',linewidth=0.75, label='$m = %.3f \pm %.3f$'%(m, np.sqrt(cov_m)))
    #axs[0].fill_between(ds_fit, model(ds_fit,xi) -err_model(ds_fit,xi,np.sqrt(cov_xi)), model(ds_fit,xi) +err_model(ds_fit,xi,np.sqrt(cov_xi)), alpha=0.2, color='red')
    axs[0].fill_between(ds, 2, y2=-1,where= mask, alpha=0.2, color='red',label = 't mins')
    axs[0].fill_between(ds, 2, y2=-1,where= mask2, alpha=0.2, color='blue', label = 't maxs')
    lows = [lowers for i in range(len(uppers))]
    Probs = Probs.reshape(len(uppers),len(lowers))
    a0s = a0s.reshape(len(uppers),len(lowers))
    sigmassqs = sigmassqs.reshape(len(uppers),len(lowers))
    
    for i in range(len(uppers)):
        axs[1].plot(lows[i],Probs[i],'o',label = 't max = ' + str(uppers[i]))
        axs[2].errorbar(lows[i],a0s[i],yerr=np.sqrt(sigmassqs[i]),fmt = '.',label = 't max = ' + str(uppers[i]))

    

    axs[0].set_ylim(m-m/10,m+m/10)
    plt.xlabel(r'wall separation $d$ [$a$]')
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    axs[0].set_ylabel('wall wall correlation $C_{ww}(d)$')
    axs[1].set_ylabel('Probability')

    axs[0].legend(prop={'size':12}, frameon=True, loc='upper right')
    axs[1].legend(prop={'size':12}, frameon=True, loc='upper right')

   

    if accel == False:

        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+' manual Accel.svg'
    plt.savefig(file_name)
    plt.show()

def plot_iats_susceptibility(betas,Ns,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Susceptibility iat'
    iats_acc, err_iats_acc = np.zeros(len(betas)), np.zeros(len(betas))
    iats_no_acc, err_iats_no_acc = np.zeros(len(betas)), np.zeros(len(betas))
    for i,beta in enumerate(betas):
        file_name_no_acc= "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(Ns[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
        file_name_acc = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(Ns[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
        iat_acc, err_iat_acc = np.load(file_name_acc)
        iat_no_acc, err_iat_no_acc = np.load(file_name_no_acc)
        iats_acc[i] = iat_acc
        iats_no_acc[i] = iat_no_acc
        err_iats_acc[i] = err_iat_acc
        err_iats_no_acc[i] = err_iat_no_acc
    mass, err_mass = np.zeros(len(betas)), np.zeros(len(betas))
    for i,beta in enumerate(betas):
        observable_name = "mass state 1"
        file_name ="ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(Ns[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "   + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
        
        x,err_mass[i] = np.load(file_name)
        mass[i] = x
#N_tau,_ = Chiral_run.load_calibration(betas1[0],N[0],SU,order,N_order,accel=True)
    x = 1/mass
    
    err_x = err_mass/mass**2
    
    plt.errorbar(x,iats_acc,yerr=err_iats_acc,fmt='.',capsize=2,color='black')
    plt.yscale('log')
    plt.xscale('log')
    plt.errorbar(x,iats_no_acc,yerr=err_iats_no_acc,fmt='.',capsize=2,color='red')
    plt.yscale('log')
    plt.xscale('log')
    log_a = np.log(x)
    log_iat_acc = np.log(iats_acc)
    err_log = err_iats_acc/iats_acc
    def linear(x, z, b):
        return z*x + b
   
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat_acc,sigma=err_log, absolute_sigma=True)
    fitted_acc = [a**popt[0]*np.exp(popt[1]) for a in x]
    chis_sq = np.sum(((log_iat_acc - linear(log_a,*popt))/err_log)**2)
    print(chis_sq/(len(x)-1))
    plt.plot(x,fitted_acc,color = 'black',linestyle='dashed',label = "Fit for accelration: $z = $" + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)))
    log_iat_no_acc = np.log(iats_no_acc)
    err_log_no_acc = err_iats_no_acc/err_iats_no_acc
    def linear(x, z, b):
        return z*x + b
   
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat_no_acc,sigma=err_log_no_acc, absolute_sigma=True)
    fitted_acc = [a**popt[0]*np.exp(popt[1]) for a in x]
    chis_sq = np.sum(((log_iat_no_acc - linear(log_a,*popt))/err_log_no_acc)**2)
    print(chis_sq/(len(x)-1))
    plt.plot(x,fitted_acc,color = 'red',linestyle='dashed',label = "Fit for No accelration: $z = $" + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)))
    
    plt.legend()
    plt.savefig("ChiralResults/Processed/Plots/Susceptibility iat" + " SU = " +str(SU)+".svg")   

    plt.show()

def plot_SU2_massoverlambda(betas,Ns,order,N_order,N_measure,N_thermal):
    SU = 2
    
    corrlens = np.zeros(len(betas))
    err_corrlens = np.zeros(len(betas))
    for i,beta in enumerate(betas):
        observable_name = "correlation length"
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(Ns[i])  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "   + str(N_measure)+" N Thermal = "  + str(N_thermal)+' Accel.npy'
        corrlens[i],err_corrlens[i] = np.load(file_name)
    print(corrlens)
    N = 2
    b0 = N / (8*np.pi)
    b1 = N**2 / (128*np.pi**2)
    G1 = 0.04616363
    b2 = 1/(2*np.pi)**3 * (N**3)/128 * ( 1 + np.pi*(N**2 - 2)/(2*N**2) - np.pi**2*((2*N**4-13*N**2+18)/(6*N**4) + 4*G1) ) 
    pre_factor = (2*np.pi*betas)**(1/2) * np.exp(-2*np.pi*betas)
    from scipy.integrate import quad
    def integrand(x):
        '''integrand in expression for the renormalisation scale using the beta function at 3 loop accuracy.
        
        Parameters
        ----------
        x: float
            value of the coupling constant squared i.e. x=g^2

        Returns
        -------
        inte: float
            the integrand
        '''
        beta_3l = -b0*x**2 - b1*x**3 - b2*x**4  
        inte = 1/beta_3l + 1/(b0*x**2) - b1/(b0**2*x)
        return inte
    
    F = np.zeros_like(betas)

    for i,beta in enumerate(betas):
        res, err = quad(integrand, 0, 4/(N*beta))
        F[i] = pre_factor[i] * np.exp(-res)

    mass_lambda_int = 1/corrlens * 1/F
    mass_lambda_int_err = mass_lambda_int / corrlens * err_corrlens
    
    #mass_lambda = 1/corrlens *  np.exp(2*np.pi*betas) / np.sqrt(2*np.pi*betas)
    #mass_lambda_err = mass_lambda / corrlens * err_corrlens
    cts_prediction = 32 * np.exp(np.pi/4) / np.sqrt(np.pi*np.e)
    fig = plt.figure(figsize=(8,6))
    plt.errorbar(betas, mass_lambda_int, yerr=mass_lambda_int_err, fmt='.', capsize=2,color = 'k')
    plt.hlines(cts_prediction, betas[0], betas[-1], linestyles='--', color='k')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$M / \Lambda_{L,2l}$')
    plt.savefig("ChiralResults/Processed/Plots/mass_lambda_int" + " SU = " +str(SU)+".svg")
    plt.show()
def plot_CV(betas, model_params, accel = False):    
    """
    Function used to plot generic data of an observable
    against beta for different parameters such as
    the lenght of the lattice or the group used
    """
    observable_name = "C_V"
    symbol = "$C_V$"

    result = [0 for i in range(len(betas))]

    err = [0 for j in range(len(betas))]

    for i, beta in enumerate(betas):
        if accel == False:

            file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                        observable_name + " beta = " + str(beta) + " N = " + \
                        str(model_params["n_length"]) + \
                        " SU = " + str(model_params["su_parameter"]) + " Order = " \
                        + str(model_params["order"]) + " N Order = " \
                        + str(model_params["n_order"]) + " N measurements = " + \
                        str(model_params["n_measure"]) + \
                        " N Thermal = " + str(model_params["n_thermal"]) + '.npy'
        else:
            file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                        observable_name + " beta = " + str(beta) + " N = " + \
                        str(model_params["n_length"]) + \
                        " SU = " + str(model_params["su_parameter"]) + " Order = " \
                        + str(model_params["order"]) + " N Order = " \
                        + str(model_params["n_order"]) + " N measurements = " + \
                        str(model_params["n_measure"]) + \
                        " N Thermal = " + str(model_params["n_thermal"]) + ' Accel.npy'

        result[i], err[i] = np.load(file_name)

    axis = plt.subplot(111)

    plt.xlabel('$\\beta$')

    plt.ylabel(symbol)

    axis.spines[['right', 'top']].set_visible(False)

    plt.xlim(0, max(betas) + 0.1)

    axis.errorbar(x=betas, y=result, yerr=err, fmt="xk", label='Data with Taylor')
    axis.plot(betas,result,'k',linestyle='dashed',linewidth = 0.5)
    plt.show()
    plt.savefig("ChiralResults/Processed/Plots/Gen Obs"+".svg")

def plot_CV_mult(betas, model_params, SUs,N_measures,N_thermals,orders,N_orders, accel = False):    
    """
    Function used to plot generic data of an observable
    against beta for different parameters such as
    the lenght of the lattice or the group used
    """
    observable_name = "C_V"
    symbol = "$C_V$"
    axis = plt.subplot(111)
    for l,SU in enumerate(SUs):
        print(SU)

        result = [0 for i in range(len(betas[l]))]

        err = [0 for j in range(len(betas[l]))]

        for i, beta in enumerate(betas[l]):
            if accel == False:

                file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                            observable_name + " beta = " + str(beta) + " N = " + \
                            str(model_params["n_length"]) + \
                            " SU = " + str(SUs) + " Order = " \
                            + str(orders[l]) + " N Order = " \
                            + str(N_orders[l]) + " N measurements = " + \
                            str(N_measures[l]) + \
                            " N Thermal = " + str(N_thermals[l]) + '.npy'
            else:
                file_name = "ChiralResults/Processed/" + observable_name + "/" + \
                            observable_name + " beta = " + str(beta) + " N = " + \
                            str(model_params["n_length"]) + \
                            " SU = " + str(SU) + " Order = " \
                            + str(orders[l]) + " N Order = " \
                            + str(N_orders[l]) + " N measurements = " + \
                            str(N_measures[l]) + \
                            " N Thermal = " + str(N_thermals[l]) + ' Accel.npy'

            result[i], err[i] = np.load(file_name)

        

        plt.xlabel('$\\beta$')

        plt.ylabel(symbol)

        axis.spines[['right', 'top']].set_visible(False)

        #plt.xlim(0, max(betas) + 0.1)
        print(result)
        axis.errorbar(x=np.array(betas[l]), y=np.array(result).copy(), yerr=err, fmt="xk", label='Data with Taylor')
        axis.plot(betas[l],result,'k',linestyle='dashed',linewidth = 0.5)
    
    plt.savefig("ChiralResults/Processed/Plots/Gen Obs"+".svg")
    plt.show()
    