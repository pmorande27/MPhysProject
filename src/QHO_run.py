from QSHOHMC import Lattice
import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats
from scipy.optimize import curve_fit


def analytical(a,N,mu,m):

    R = 1+ a**2*mu**2/(2*m)-a*mu/np.sqrt(m)*(1+a**2*mu**2/(4*m))**0.5

    return 1/(2*mu*(m+a**2*mu**2/4)**0.5)*(1+R**N)/(1-R**N)
def main():
    
   
    
    pos_a =np.logspace(0,-2,20)
    #print(pos_a)
    #pos_a = [0.01930698]
    #pos_a = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    
    #pos_a = [0.6]
    #pos_a.reverse()
    #pos_a =  [0.1]
    N = 1000
    mass = 1
    w = 1
    N_measurements = 10**5

    #lat = Lattice(100,0.1,1000,10**5,1,1,N_tau,d_tau,1,1,True)
    #lat.generate_measurements(Lattice.measure_position)
    #print(np.mean(np.exp(-lat.DH)))
    
    [calibration(a,N,mass,w,2,True) for a in pos_a]

    #[calibration(a,N,mass,w,3,True) for a in pos_a]
    #[measure_two_point_function_a(a,N,N_measurements,acceleration=False) for a in pos_a]
    
    #[measure_two_point_function_a(a,N,N_measurements,acceleration=True) for a in pos_a]
    #measure_DH(0.5,N,10**5,False)
    #measure_sq_a(pos_a,N,N_measurements,True)
    #plot_sq_a(pos_a,N,N_measurements,True)

    #plot_DH(0.5,N,10**5,False) 
    #measure_DH(0.5,N,10**5,True)
    #plot_DH(0.5,N,10**5,True)
    #plot_position_sq_accel_and_no_accel(pos_a,N,N_measurements)
    #plot_sq_a(pos_a,N,N_measurements,True)
    #measure_sq_a(pos_a,N,N_measurements,True)
    #[obtain_model_2(a,N,N_measurements,True)for a in pos_a]
    #[measure_config_a(a,N,N_measurements,False) for a in pos_a]
    #[measure_iat(a,N,N_measurements,True) for a in pos_a]
    #plot_iat_accel_and_no_accel(pos_a,N,N_measurements)

    #plot_models_2(pos_a,N,N_measurements,acceleration=True)
    #plot_iat_estimation(0.2,100,False)
def plot_DH(a, N, N_measurements, acceleration):
    file_name = "QHOResults/IAT/Measure DH "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = "+str(acceleration) +".npy"
    DH = np.load(file_name)
    average = np.average(np.exp(-DH))
    ax = plt.subplot(111)

    ax.plot(np.exp(-DH),'.k',label = 'Data')
    ax.plot([0,len(DH)],[average,average],label = 'Average value')
    plt.xlim(500,1800)
    plt.legend()
    print(average)
    print((np.var(np.exp(-DH),ddof=1))/len(DH)**0.5)

    ax.spines[['right', 'top']].set_visible(False)
    plt.xlabel("Run")
    plt.ylabel('$ \exp\{-\Delta H\} $')
    plt.savefig('QHOResults/Plots/ExpDH_accel_' +str(acceleration)+ '.svg')


    plt.show()

def measure_DH(a, N, N_measurements, acceleration):
    N_tau,d_tau = load_calibration(a,N,1,1,acceleration)
    lat = Lattice(N,a,int(N_measurements/10),N_measurements,1,1,N_tau,d_tau,1,1,acceleration)
    values = lat.generate_measurements(observable=Lattice.measure_position)
    Hs = lat.DH
    file_name = "QHOResults/IAT/Measure DH "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = "+str(acceleration) 
    np.save(file_name,Hs)

def plot_iat_estimation(a, N,acceleration):
    N_tau, d_tau = load_calibration(a,N,1,1,acceleration)
    Ns = [int(i) for i in np.logspace(6,3,20)]

    IATs = np.zeros(len(Ns))
    IATsErr = np.zeros(len(Ns))
    for i,N_measurement in enumerate(Ns):
        lat = Lattice(N, a, int(N_measurement/10), N_measurement,1,1,N_tau,d_tau,1,1,acceleration)
        values = np.array(lat.generate_measurements(Lattice.measure_position))
        IATs[i], IATsErr[i] = Stats.autocorrelator(values)
    plt.plot(Ns,IATs)    
    plt.xscale('log')
    plt.show()

        
def plot_two_point_function_accel_and_no_accel(a, N,N_measurements):
    file_name_acc = "QHOResults/Two Point Function/Measure Two Point Function "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = "+str(True) +".npy"
    file_name_noacc = "QHOResults/Two Point Function/Measure Two Point Function "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = "+str(False) +".npy"
    vals_acc = np.load(file_name_acc)

    vals_noacc = np.load(file_name_noacc)
    Stat_noacc = [Stats(b) for b in vals_noacc]
    
    estimates_noacc = [S.estimate() for S in Stat_noacc]
    
    value_estimates_noacc = [estimates_noacc[i][0] for i in range(N+1)]
    
    error_estimates_noacc = [estimates_noacc[i][1] for i in range(N+1)]
    
    normalized_data_noacc, normalized_error_noacc = value_estimates_noacc/value_estimates_noacc[0], error_estimates_noacc/value_estimates_noacc[0] 
    
    Stat_acc = [Stats(b) for b in vals_acc]
    
    estimates_acc = [S.estimate() for S in Stat_acc]
    
    value_estimates_acc = [estimates_acc[i][0] for i in range(N+1)]
    
    error_estimates_acc = [estimates_acc[i][1] for i in range(N+1)]
    
    normalized_data_acc, normalized_error_acc = value_estimates_acc/value_estimates_acc[0], error_estimates_acc/value_estimates_acc[0] 
    
    ax = plt.subplot(111)

    ax.errorbar([i for i in range(N+1)],normalized_data_acc,yerr=normalized_error_acc,fmt='.k', label = 'Simulation Data Accelerated')
    ax.errorbar([i for i in range(N+1)],normalized_data_noacc,yerr=normalized_error_noacc,fmt='.k', label = 'Simulation Data Not Accelerated',color = 'red')

    
    ax.spines[['right', 'top']].set_visible(False)

    plt.show()

    

def plot_iat_accel_and_no_accel(pos_a, N, N_measurements):
    iats_acc = np.zeros(len(pos_a))
    iats_acc_err = np.zeros(len(pos_a))
    iats_noacc = np.zeros(len(pos_a))
    iats_noacc_err = np.zeros(len(pos_a))
    ax = plt.subplot(111)


    
    for i in range(len(pos_a)):
        file_name_acc = "QHOResults/IAT/Measure IAT "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i])+ " Accel = "+str(True) +".npy"
        file_name_noacc = "QHOResults/IAT/Measure IAT "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i])+ " Accel = "+str(False) +".npy"

        value_acc = np.load(file_name_acc)
        iats_acc[i] = value_acc[0]
        iats_acc_err[i] = value_acc[1]
        value_noacc = np.load(file_name_noacc)
        iats_noacc[i] = value_noacc[0]
        iats_noacc_err[i] = value_noacc[1]
    log_a = np.log(pos_a)
    log_iat = np.log(iats_noacc)
    err_log = iats_noacc_err/iats_noacc
    def linear(x, z, b):
        return z*x + b
    popt, pcov = curve_fit(linear, xdata=log_a,ydata=log_iat,sigma=err_log, absolute_sigma=True)
    fitted_noacc = [a**popt[0]*np.exp(popt[1]) for a in pos_a]
    ax.errorbar(x=pos_a, y = iats_acc, yerr = iats_acc_err, fmt = '.k', color = 'red',label = "Acceleration")
    ax.errorbar(x=pos_a, y = iats_noacc, yerr = iats_noacc_err, fmt = '.k',label = "No Acceleration")
    log_iat_2 = np.log(iats_acc)
    err_log_2 = iats_acc_err/iats_acc
    popt2, pcov2 = curve_fit(linear, xdata=log_a,ydata=log_iat_2,sigma=err_log_2, absolute_sigma=True)
    fitted_acc = [a**popt2[0]*np.exp(popt2[1]) for a in pos_a]
    ax.plot(pos_a,fitted_noacc,label = 'Fit for no accelration: $z = $' + str(round(popt[0],3))+' $\pm $'+ str(round(np.sqrt(pcov[0][0]),3)))

    ax.plot(pos_a,fitted_acc,label = 'Fit for accelration: $z = $' + str(round(popt2[0],3))+' $\pm $'+ str(round(np.sqrt(pcov2[0][0]),3)))


    print(popt[0], np.sqrt(pcov[0][0]))
    print(popt2[0], np.sqrt(pcov2[0][0]))
    plt.xscale('log')
    plt.yscale('log')
    ax.spines[['right', 'top']].set_visible(False)
    plt.ylabel('$ \\tau_{INT}$')
    plt.xlabel('a')
    plt.legend()

    plt.savefig('QHOResults/Plots/Iat_Accel_no_AccelModified_HMC_noaccel_1' + '.svg')

    plt.show()

def plot_position_sq_accel_and_no_accel(pos_a,N, N_measurements):
    result_acc = np.zeros(len(pos_a))

    error_acc = np.zeros(len(pos_a))

    error_noacc = np.zeros(len(pos_a))
    
    result_noacc = np.zeros(len(pos_a))

    for i in range(len(pos_a)):

        file_name_acc = "QHOResults/Sq Position/Measure Sq Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " +str(True)+".npy"
        file_name_noacc = "QHOResults/Sq Position/Measure Sq Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " +str(False)+".npy"


        values_acc = np.load(file_name_acc)

        result_acc[i],error_acc[i] = Stats(values_acc).estimate()

        values_noacc = np.load(file_name_noacc)
        result_noacc[i], error_noacc[i] = Stats(values_noacc).estimate()
    real = [analytical(a,N,1,1) for a in np.linspace(0.1,1,100)]

    ax = plt.subplot(111)
    
    ax.errorbar(x=pos_a,y=result_acc,yerr= error_acc,fmt=".k",label = 'Accelerated Data')
    ax.errorbar(x=pos_a,y=result_noacc,yerr= error_noacc,fmt=".k",label = 'Not Accelerated Data',color = "red")

    ax.plot(np.linspace(0.1,1,100),real,"black",linestyle=(0, (1, 1)),label = "Discrete Theory")

    ax.plot([0,1.5],[0.5,0.5],color = "black", label = "Continuum Limit",linestyle="--")

    plt.xlim(0,1.1)
    
    plt.legend()
    
    plt.xlabel('$\epsilon$')
    
    plt.ylabel('$\langle x^2\\rangle$')
    
    ax.spines[['right', 'top']].set_visible(False)
    
    plt.savefig('QHOResults/Plots/Position_sq_Accel_no_Accel' + '.svg')

    plt.show()

def plot_position_accel_and_no_accel(pos_a, N, N_measurements):
    result_acc = np.zeros(len(pos_a))

    error_acc = np.zeros(len(pos_a))

    error_noacc = np.zeros(len(pos_a))
    
    result_noacc = np.zeros(len(pos_a))

    for i in range(len(pos_a)):

        file_name_acc = "QHOResults/Position/Measure Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " +str(True)+".npy"
        file_name_noacc = "QHOResults/Position/Measure Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " +str(False)+".npy"


        values_acc = np.load(file_name_acc)

        result_acc[i],error_acc[i] = Stats(values_acc).estimate()

        values_noacc = np.load(file_name_noacc)
        result_noacc[i], error_noacc[i] = Stats(values_noacc).estimate()

    ax = plt.subplot(111)
    
    ax.errorbar(x=pos_a,y=result_acc,yerr= error_acc,fmt=".k",label = 'Accelerated Data')
    ax.errorbar(x=pos_a,y=result_noacc,yerr= error_noacc,fmt=".k",label = 'Not Accelerated Data',color = 'red')

    ax.plot([0,1.5],[0,0],color = "black", label = "Continuum Limit",linestyle="--")

    plt.xlim(0,1.1)
    
    plt.legend()
    
    plt.xlabel('$\epsilon$')
    
    plt.ylabel('$\langle x\\rangle$')
    
    ax.spines[['right', 'top']].set_visible(False)

    plt.savefig('QHOResults/Plots/Position_Accel_no_Accel' + '.svg')

    plt.show()






    
def measure_iat(a,N, N_measurements, acceleration):


    mass = 1
    w = 1
    N_tau , d_tau = load_calibration(a, N, mass, w,acceleration)
    N_thermal = int(N_measurements/10)
    N_sweeps = 1
    Delta = 1

    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=N_tau,d_tau=d_tau,mass=mass,w=w,acceleration=acceleration)
    values = lat.generate_measurements(Lattice.measure_position)
    iat,error = Stats.autocorrelator(np.array(values))
    file_name = "QHOResults/IAT/Measure IAT "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = "+str(acceleration) 
    np.save(file_name,np.array([iat,error]))

def plot_iat(pos_a, N, N_measurements, acceleration):
    iats = np.zeros(len(pos_a))
    
    for i in range(len(pos_a)):
        file_name = "QHOResults/IAT/Measure IAT "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i])+ " Accel = "+str(acceleration) +".npy"
        value = np.load(file_name)
        iats[i] = value[0]
    plt.plot(pos_a,iats,'.k')
    plt.xscale('log')
    plt.yscale('log')


    
    
    
    
def measure_sq_a(pos_a,N,N_measurements,acceleration= False):

    for i in range(len(pos_a)):    
        
    
        a = pos_a[i]
    
    
        N_sweeps = 1
    
        N_thermal = int(N_measurements/10)
    
        Delta = 1

        w = 1

        mass = 1

        N_tau , d_tau = load_calibration(a, N, mass, w,acceleration=acceleration)

        lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=N_tau,d_tau=d_tau,mass=mass,w=w,acceleration=acceleration)
    
        values = lat.generate_measurements(Lattice.measure_sq_position)

        file_name = "QHOResults/Sq Position/Measure Sq Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a) + " Accel = " +str(acceleration)
        print(Stats(values).estimate())
        np.save(file_name,values)

def plot_sq_a(pos_a, N, N_measurements,acceleration=False):

    result = [0 for i in range(len(pos_a))]

    error = [0 for i in range(len(pos_a))]

    for i in range(len(pos_a)):

        file_name = "QHOResults/Sq Position/Measure Sq Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " +str(acceleration)+".npy"

        values = np.load(file_name)
        

        result[i],error[i] = Stats(values).estimate()
        print(pos_a[i],result[i],error[i])

    real = [analytical(a,N,1,1) for a in np.linspace(min(pos_a),1,1000)]

    ax = plt.subplot(111)
    
    ax.errorbar(x=pos_a,y=result,yerr= error,fmt=".k")
    
    ax.plot(np.linspace(min(pos_a),1,1000),real,"black",linestyle=(0, (1, 1)),label = "Discrete Theory")

    ax.plot([0,1.5],[0.5,0.5],color = "black", label = "Continuum Limit",linestyle="--")

    plt.xscale('log')
    plt.xlim(0,1.1)
    
    plt.legend()
    
    plt.xlabel('$\epsilon$')
    
    plt.ylabel('$\langle x^2\\rangle$')
    
    ax.spines[['right', 'top']].set_visible(False)

    plt.savefig('QHOResults/Plots/QHO_sq_diff_a_Accel_'+str(acceleration)+  '.svg')
    
    plt.show()

    
def measure_two_point_function_a(a, N, N_measurements, acceleration = False):
    
    N_sweeps = 1

    N_thermal =int(N_measurements/10)

    Delta = 1

    w = 1
    
    m = 1
    
    N_tau,d_tau = load_calibration(a,N,m,w,acceleration)


    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=N_tau,d_tau=d_tau, mass=m, w=w, acceleration=acceleration)
    
    values = lat.generate_measurements(Lattice.measure_greens)

    vals = [[values[i][n] for i in range(N_measurements)] for n in range(N+1)]

    file_name = "QHOResults/Two Point Function/Measure Two Point Function " + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = " +str(acceleration)
    
    np.save(file_name,vals)

def plot_two_point_function_a(a, N, N_measurements,acceleration = False):

    file_name = "QHOResults/Two Point Function/Measure Two Point Function " + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = " +str(acceleration)+".npy"
    
    vals = np.load(file_name)
    
    Stat = [Stats(b) for b in vals]
    
    estimates = [S.estimate() for S in Stat]
    
    value_estimates = [estimates[i][0] for i in range(N+1)]
    
    error_estimates = [estimates[i][1] for i in range(N+1)]
    
    normalized_data, normalized_error = value_estimates/value_estimates[0], error_estimates/value_estimates[0] 
    
    ax = plt.subplot(111)

    ax.errorbar([i for i in range(N+1)],normalized_data,yerr=normalized_error,fmt='.k', label = 'Simulation Data')
    
    ax.spines[['right', 'top']].set_visible(False)

    plt.show()

def obtain_model_a(a, N, N_measurements, acceleration = False, plot = False):
    
    file_name = "QHOResults/Two Point Function/Measure Two Point Function " + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a) + " Accel = "+ str(acceleration) + ".npy"
    
    vals = np.load(file_name)
    
    Stat = [Stats(b) for b in vals]
    
    estimates = [S.estimate() for S in Stat]
    
    value_estimates = [estimates[i][0] for i in range(N+1)]
    
    error_estimates = [estimates[i][1] for i in range(N+1)]
    
    def model(t,dE):
    
        return (np.cosh(dE*(t-N/2)) ) / (np.cosh(dE*N/2) )
    
    normalized_data, normalized_error = value_estimates[:-1]/value_estimates[0], error_estimates[:-1]/value_estimates[0] 
    
    sep = np.arange(N)
    
    upper = 20
    
    mask = sep <= upper
    
    sep_fit = sep[mask]
    
    popt, pcov = curve_fit(model, sep_fit, normalized_data[mask], sigma=normalized_error[mask], absolute_sigma=True, bounds=(0,np.inf))

    if plot:
        ax = plt.subplot(111)

        ax.errorbar([i for i in range(N)],normalized_data,yerr=normalized_error,fmt='.k', label = 'Simulation Data')
        
        xs = np.linspace(0,N)
        
        ys = [model(x,popt[0]) for x in xs]
        
        ax.plot(xs,ys,"black",linestyle =(0, (1, 1)),label = "Best Fit" )
        
        plt.xlabel('t')
        
        plt.ylabel("$\langle x_{t_0} x_{t_0+t} \\rangle$")
        
        plt.legend()
        
        ax.spines[['right', 'top']].set_visible(False)
        
        plt.savefig('QHOResults/QHO_twopointModel a = '+str(a)+'.svg')

        plt.show()
    return popt[0]/a ,np.sqrt(pcov[0][0]) / a

def measure_pos_a(pos_a,N,N_measurements, acceleration):
    for i in range(len(pos_a)):    
    
        a = pos_a[i]

        N_sweeps = 1
    
        N_thermal = int(N_measurements/10)
    
        Delta = 1
        
        mass = 1
        
        w = 1
        
        N_tau , d_tau = load_calibration(a,N,mass,w, acceleration)
    
        lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=N_tau,d_tau=d_tau,mass=mass,w=w,acceleration=acceleration)
    
        values = lat.generate_measurements(Lattice.measure_position)
        
        file_name = "QHOResults/Position/Measure Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a) + " Accel = " + str(acceleration)
        
        np.save(file_name,values)

def plot_pos_a(pos_a, N, N_measurements, acceleration):
    
    plt.xlabel('a')

    result = [0 for i in range(len(pos_a))]
    
    error = [0 for i in range(len(pos_a))]
    for i in range(len(pos_a)):
        
        file_name = "QHOResults/Position/Measure Position "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i]) + " Accel = " + str(acceleration) +".npy"
        
        values = np.load(file_name)
        
        result[i],error[i] = Stats(values).estimate()
        
    ax = plt.subplot(111)
    
    ax.errorbar(x=pos_a,y=result,yerr=error,fmt=".k")
    
    ax.plot(np.linspace(0,1.1,100),[0 for i in range(100)],"black",linestyle=(0, (1, 1)),label = "Theory")
    
    
    plt.legend()
    plt.xscale('log')
    plt.xlim(0,1.1)

    
    plt.xlabel('$\epsilon$')
    
    plt.ylabel('$\langle x\\rangle$')
    
    ax.spines[['right', 'top']].set_visible(False)

    plt.savefig('QHOResults/Plots/QHO_pos_diff_a_accel_Accel_' +str(acceleration) + '.svg')
    
    plt.show()

def measure_config_a(a,N,N_measurements, acceleration = False):
    N_sweeps = 1
    
    N_thermal =int(N_measurements/10)
    
    Delta = 1
    
    mass = 1
    
    w = 1
    
    N_tau, d_tau = load_calibration(a,N,mass,w,acceleration)
    

    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau= N_tau,d_tau = d_tau, mass=mass,w=w, acceleration=acceleration)
    
    phi = lat.Get_configurations()
    
    file_name = "QHOResults/Configurations/Measure Configurations "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+" Accel = " + str(acceleration)
    
    np.save(file_name,phi)

def plot_wvfn(a,N,N_measurements,acceleration=False):
    
    values =     np.load("QHOResults/Configurations/Measure Configurations "  + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+" Accel = " +str(acceleration)+".npy")
    
    values =values.flatten()
    
    ax = plt.subplot(111)

    ax.hist(values,bins=500,color='grey',density=True,alpha = 0.4)
    
    xs = np.linspace(-4,4,1000)
    
    w = (1+a**2/4)**0.5
    
    R = 1+ a**2/(2)-a*(1+a**2/(4))**0.5
    
    w2 = (1+a**2/4)**0.5 *(1-R**N)/(1+R**N)
    
    ys = [1/np.sqrt(np.pi)*np.exp(-x**2) for x in xs]
    
    ys_2 =  [w**0.5/np.sqrt(np.pi)*np.exp(-w*x**2) for x in xs]
    
    ys_3 =[w2**0.5/np.sqrt(np.pi)*np.exp(-w2*x**2) for x in xs]
    
    ax.spines[['right', 'top']].set_visible(False)
    
    ax.plot(xs,ys,'black',linestyle = '--',label = "Continuum Limit")
    
    ax.plot(xs,ys_2,'black',linestyle = ':',label = "Discrete Theory")
    
    ax.plot(xs,ys_3,'red',linestyle = ':',label = "Discrete Theory 2")

    
    plt.xlabel('$x$')
    
    plt.ylabel("$|\psi(x)|^2$")
    
    plt.legend()
    
    plt.savefig('QHOResults/Plots/QHO_wvfn a=1.svg')

    plt.legend()



    plt.show()
def plot_models(pos_a,N,N_measurements,acceleration = False):

    values = [0 for i in range(len(pos_a))]

    errors = [0 for i in range(len(pos_a))]

    for i in range(len(pos_a)):
       
        values[i],errors[i] = obtain_model_a(pos_a[i],N,N_measurements,acceleration=acceleration,plot=False)
    
    ax =plt.subplot(111)
    
    ax.errorbar(pos_a,values,yerr=errors ,fmt='.k', label = 'Simulation Data')
    
    plt.xlabel('a')
    
    plt.ylabel('$E_1-E_0$')
    
    mu= 1
    
    m =1
    
    j = 1
    def A(a):
        return 1*np.sqrt(1+1/4*(a)**2)
    def R(a):
        return np.sqrt(1+(A(a)*a)**2)- a*A(a)
    
    j =5
    
    xs = np.linspace(0.1,1,100)
    
    ys2 = [-1/a * np.log((R(a)**(j+1)-R(a)**(N-j-1))/(R(a)**j-R(a)**(N-j)))for a in xs]
    
    ax.plot(xs,ys2,color = 'black', linestyle = ':',label = 'Discrete Theory')
    
    ax.spines[['right', 'top']].set_visible(False)
    
    ax.plot([0,1.2],[1,1],color = 'black', linestyle = '--',label = 'Continuum Limit')
    
    plt.xlim(0,1.1)
    
    plt.legend()

    plt.savefig('QHOResults/Plots/QHO_DeltaE.svg')
    
    
    plt.show()

    #plt.plot(xs,ys)
    





def obtain_model_2(a, N, N_measurements, acceleration=False):
    N_sweeps = 1

    N_thermal =int(N_measurements/10)

    Delta = 1

    w = 1
    
    m = 1
    N_tau, d_tau = load_calibration(a,N,m,w,acceleration)

    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=N_tau,d_tau=d_tau, mass=m, w=w, acceleration=acceleration)
    values = lat.generate_measurements(Lattice.measure_greens)

    vals = [[values[i][n] for i in range(N_measurements)] for n in range(N+1)]
    file_name = "QHOResults/Ediff/Delta E " + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(a)+ " Accel = " + str(acceleration)+".npy"
    
    Stat = [Stats(b) for b in vals]
    estimates = [S.estimate() for S in Stat]
    value_estimates = [estimates[i][0] for i in range(N+1)]
    error_estimates = [estimates[i][1] for i in range(N+1)]
    def model(t,dE):
        return (np.cosh(dE*(t-N/2)) - 1) / (np.cosh(dE*N/2) - 1)
    normalized_data, normalized_error = value_estimates[:-1]/value_estimates[0], error_estimates[:-1]/value_estimates[0] 
    sep = np.arange(N)
    upper = 20
    mask = sep <= upper
    sep_fit = sep[mask]
    popt, pcov = curve_fit(model, sep_fit, normalized_data[mask], sigma=normalized_error[mask], absolute_sigma=True, bounds=(0,np.inf))
    np.save( file_name,[popt[0]/a ,np.sqrt(pcov[0][0]) / a])

def plot_models_2(pos_a, N, N_measurements, acceleration = False):
    values = [0 for i in range(len(pos_a))]
    errors = [0 for i in range(len(pos_a))]
    for i in range(len(pos_a)):
        file_name = "QHOResults/Ediff/Delta E " + "N = "+str(N) + " N_measure = " + str(N_measurements) + " a = " + str(pos_a[i])+ " Accel = " + str(acceleration)+".npy"
        val = np.load(file_name)
        values[i] = val[0]
        errors[i] = val[1]
    ax =plt.subplot(111)
    ax.errorbar(pos_a,values,yerr=errors ,fmt='.k', label = 'Simulation Data')
    plt.xlabel('$\epsilon$')
    plt.ylabel('$E_1-E_0$')
    
    mu= 1
    m =1
    
    j = 1
    def A(a):
        return 1*np.sqrt(1+1/4*(a)**2)
    def R(a):
        return np.sqrt(1+(A(a)*a)**2)- a*A(a)
    j =5
    xs = np.linspace(0.1,1,100)
    ys2 = [-1/a * np.log((R(a)**(j+1)-R(a)**(N-j-1))/(R(a)**j-R(a)**(N-j)))for a in xs]
    ax.plot(xs,ys2,color = 'black', linestyle = ':',label = 'Discrete Theory')
    ax.spines[['right', 'top']].set_visible(False)
    


    ax.plot([0,1.2],[1,1],color = 'black', linestyle = '--',label = 'Continuum Limit')
    plt.xlim(0,1.2)
    plt.legend()

    plt.savefig('QHOResults/Plots/QHO_DeltaE 2.svg')
    plt.show()

def calibration(a, N, mass, w, N_tau_guess = 2, acceleration = False):
    Delta = 1
    N_tau = N_tau_guess
    print('Calibration with a = ' + str(a) + " N = " +str(N)+ " and Accel = " + str(acceleration))

    while True:
        d_tau = 1/N_tau
        calibration_runs = 10**3

        lat = Lattice(N = N, a = a, N_thermal = 0, N_measurement = 0, N_sweeps = 1, Delta = 1,N_tau=N_tau,d_tau=d_tau, mass=mass, w=w,acceleration=acceleration)
        lat.Calibration_Runs(calibration_runs, 1000)
        rate = lat.accepted/lat.tries
        print(rate,N_tau)
        if rate <=0.85 and rate >= 0.55:
            break
        if rate < 0.85:
            N_tau += 1
        else:
            N_tau -= 1
    print("-----------------")
    print(rate)
    file_name = "QHOParams/QHO Calibration parameters a = " + str(a) + " N = " + str(N)  + " m = " + str(mass) + " w = "+ str(w) + " Acc = " + str(acceleration)
    print(N_tau)
    np.save(file_name, [N_tau,1/N_tau])
def load_calibration(a, N, mass, w, acceleration = False):
    file_name = "QHOParams/QHO Calibration parameters a = " + str(a) + " N = " + str(N)  + " m = " + str(mass) + " w = "+ str(w) + " Acc = " + str(acceleration)+'.npy'
    parameters = np.load(file_name)
    N_tau, d_tau = int(parameters[0]), parameters[1]
    return N_tau, d_tau
main()