"""
Module used to make all hte relevants plots for the Chiral Model
"""
import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats
from scipy.optimize import curve_fit
def mass_plot(beta,N,SU,order,N_order,N_measure,N_thermal,lower_c,upper_c,lower_e,upper_e,accel = False):
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

    delta_x = 1/2 * (A*(rel_err_1 - rel_err) + B*(rel_err__1 - rel_err))
    # delta_x = A/2*(np.sqrt(rel_err_1**2 + rel_err**2)) + B/2*(np.sqrt(rel_err__1**2 + rel_err**2))
    m_eff_err = 1/np.sqrt(x**2-1) * delta_x

    cor_1 = np.roll(cor, -1) # shift to d+1
    #m_eff = - np.log(cor_1 / cor)
    #m_eff_err = np.roll(cor_err, -1)/cor_1 - cor_err/cor 
    #plt.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.', capsize=2)
    def model(x, m):
        return m
    one = ds <= upper_c
    two = ds >= lower_c

    mask = one & two
    popt, pcov = curve_fit(model, ds[mask], m_eff[mask], sigma=m_eff_err[mask], absolute_sigma=True)
    cor_len = 1/popt[0]
    cor_len_err = np.sqrt(pcov[0][0])/popt[0]**2
    ys = np.array([popt[0] for x in ds[mask]])
    r = m_eff[mask] - ys
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    reduced_chi2 = np.sum((r/m_eff_err[mask])**2) / (mask.size -1)
    axis.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.',color='black', capsize=2,label='Cosh')
    axis.plot(ds[mask], ys, c='black',linestyle='dashed',linewidth=1.0,label='$\\xi^{cosh} = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_len, cor_len_err, reduced_chi2))
    #plt.plot(ds, model(ds, *popt), c='g', label='$m = %.3f \pm %.3f$'%(popt[0], np.sqrt(pcov[0][0])))

    m_eff = - np.log(cor_1 / cor)
    m_eff_err = np.roll(cor_err, -1)/cor_1 - cor_err/cor 
    def model(x, m):
        return m
    one = ds <= upper_e
    two = ds >= lower_e
    mask = one & two
    popt, pcov = curve_fit(model, ds[mask], m_eff[mask], sigma=m_eff_err[mask], absolute_sigma=True)
    cor_len = 1/popt[0]
    cor_len_err = np.sqrt(pcov[0][0])/popt[0]**2
    ys = np.array([popt[0] for x in ds[mask]])
    r = m_eff[mask] - ys
    axis = plt.subplot(111)
    axis.spines[['right', 'top']].set_visible(False)
    reduced_chi2 = np.sum((r/m_eff_err[mask])**2) / (mask.size -1)
    axis.errorbar(ds, m_eff, yerr=m_eff_err, fmt='.', capsize=2,label='Exp')

    axis.plot(ds[mask], ys, c='red',linestyle='dashed',linewidth=1.0,label='$\\xi^{exp} = %.3f \pm %.3f$\n $\chi^2/DoF = %.3f$'%(cor_len, cor_len_err, reduced_chi2))
    plt.legend()
    plt.ylim(popt[0]-0.2, popt[0]+0.2)
    plt.title( "beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal),fontsize = 5)
    plt.ylabel('m')
    plt.xlabel('t')
    plt.xlim(0,21+1)
    if accel == False:

        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.svg'
    else:
        file_name = "ChiralResults/Processed/Plots/Mass Plot"+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.svg"
    plt.savefig(file_name)
    plt.show()
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

            file_name = "ChiralResults/Processed" + observable_name + "/" + \
                        observable_name + " beta = " + str(beta) + " N = " + \
                        str(model_params["n_length"]) + \
                        " SU = " + str(model_params["su_parameter"]) + " Order = " \
                        + str(model_params["order"]) + " N Order = " \
                        + str(model_params["n_order"]) + " N measurements = " + \
                        str(model_params["n_measure"]) + \
                        " N Thermal = " + str(model_params["n_thermal"]) + '.npy'
        else:
            file_name = "ChiralResults/Processed" + observable_name + "/" + \
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
