import numpy as np
import Stats
import Chiral_run

def process_Greens(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros((N,N))
    errors = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            values[i][j],errors[i][j],_,_ = Stats.Stats(data[i][j]).estimate()
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(values,errors)) 
    
def process_ww_corr(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_,_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"


    np.save(file_name,(values,values_err))
def process_ww_corr_state_two(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_,_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"


    np.save(file_name,(values,values_err))
def process_G_0_mom(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens 0 mom'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_,_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(values,values_err))
def process_G_diags(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'Greens diags'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_,_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(values,values_err))
def process_Greens_2(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Greens 2'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    values = np.zeros((N,N,N,N))
    errors = np.zeros((N,N,N,N))
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    values[i][j][m][n],errors[i][j][m][n],_,_ = Stats.Stats(data[i][j][m][n]).estimate()
        
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(values,errors))
def process_iats(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Susceptibility'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)

    value, error,iat, err_iat = Stats.Stats(data).estimate()
    observable_name = 'Susceptibility iat'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    print(iat,err_iat)
    np.save(file_name,(iat,err_iat))
    return iat,err_iat

def process_cost(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Susceptibility'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)  
    value, error,iat, err_iat = Stats.Stats(data).estimate()
    
    N_tau, e = Chiral_run.load_calibration(beta, N, SU, order, N_order,accel)
    cost = N_tau * iat
    cost_error = N_tau * err_iat
    observable_name = 'Cost'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(cost,cost_error))
    return iat,err_iat

def process_cost_mass(beta,N,SU,order,N_order,N_measure,N_thermal,mass,accel = True):
    observable_name = 'Susceptibility'
    
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+ " Mass = " +str(mass)+" Accel.npy"
    data = np.load(file_name)  
    value, error,iat, err_iat = Stats.Stats(data).estimate()
    
    N_tau, e = Chiral_run.load_calibration(beta, N, SU, order, N_order,accel)
    cost = N_tau * iat
    cost_error = N_tau * err_iat
    observable_name = 'Cost'
   
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Mass = " +str(mass)+" Accel.npy"
    np.save(file_name,(cost,cost_error))
    return iat,err_iat


def process_Action(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Action'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    value, error,_,_ = Stats.Stats(data).estimate()
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(value,error))

def process_Generic_observable(beta,N,SU,order,N_order,N_measure,N_thermal, observable_name,accel = False):
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    value, error,_,_ = Stats.Stats(data).estimate()
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    print(value,error)
    np.save(file_name,(value,error))

def reprocess_Greens(beta,N,SU,order,N_order,N_measure,N_thermal, N_thermal_new,accel = False):
    observable_name = 'Greens'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    data = np.load(file_name)
    N_thermal_total = N_thermal + N_thermal_new
    values = np.zeros((N,N))
    errors = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            values[i][j],errors[i][j],_,_ = Stats.Stats(np.array(data[i][j])[N_thermal_new:]).estimate()
    N_measure = N_measure - N_thermal_new
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(values,errors)) 
def reprocess_ww(beta,N,SU,order,N_order,N_measure,N_thermal, N_thermal_new, accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    data = np.load(file_name)
    N_thermal_total = N_thermal + N_thermal_new
    values = np.zeros((N+1))
    errors = np.zeros((N+1))
    for i in range(N):
        values[i],errors[i],_,_ = Stats.Stats(np.array(data[i])[N_thermal_new:]).estimate()
    values[N],errors[N] = values[0],errors[0]
    N_measure = N_measure - N_thermal_new
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
 
    np.save(file_name,(values,errors)) 

def sus_iat_process(betas,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'Susceptibility'
    iats = np.zeros(len(betas))

    for i,beta in enumerate(betas):
        if accel == False:
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
        else:
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
        data = np.load(file_name)
        val,err,iats[i],_ = Stats.Stats(data).estimate()
    return iats, betas
def sus_iat_process_N(beta,Ns,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Susceptibility'
    iats = np.zeros(len(Ns))

    for i,N in enumerate(Ns):
        if accel == False:
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
        else:
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
        data = np.load(file_name)
        val,err,iats[i],_ = Stats.Stats(data).estimate()
    return iats, Ns


def process_mass_log_jack(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def log_mass(x):
        """X = np.concatenate((x,[x[0]]))
        x_1 = np.roll(X, -1)
        mass = - np.log(x_1/X)
        N_2 = int(N/2)
        mass2 = 1/2 * (mass[:N_2+1] + abs(mass[N_2:][::-1]))
        mass2[0] = mass[0]

        return mass2"""
        N_2 = int(N/2)
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)

        return - np.log(x_1/cor)
    m_eff,m_eff_err      =  stats.jacknife_arbitrary_function_2D_2_2(N,200,log_mass)
    observable_name = 'log mass'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    
    np.save(file_name,(m_eff,m_eff_err))
    
def process_mass_cosh_jack(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def cosh_mass(x):
        N_2 = int(N/2)
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)
        return np.arccosh((x_1+x_2) / (2*cor))
    m_eff,m_eff_err      =  stats.jacknife_arbitrary_function_2D_2_2(N,100,cosh_mass)
    observable_name = 'cosh mass'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(m_eff,m_eff_err))
    

def process_ww_corr_jack(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'ww corr'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def corlen(x):
        
        
        return x
    ww_corr,ww_corr_err,      =  stats.jacknife_arbitrary_function_2D_2(N,100,corlen)
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(ww_corr,ww_corr_err))

def process_second_moment_correlation_length(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'Greens'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    second_moment,error = stats.jacknife_arbitrary_function_second_moment_correlation_length(N,100)
    observable_name = 'second moment correlation length'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    print(second_moment,error)
    np.save(file_name,(second_moment,error))
    
def process_ww_corr_state_2_jack(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def corlen(x):
        
        
        return x
    ww_corr,ww_corr_err      =  stats.jacknife_arbitrary_function_2D_2(N,100,corlen)
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(ww_corr,ww_corr_err))

def process_mass_cosh_jack_state_2(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def cosh_mass(x):
        N_2 = int(N/2)
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)
        return np.arccosh((x_1+x_2) / (2*cor))
    m_eff,m_eff_err      =  stats.jacknife_arbitrary_function_2D_2_2(N,100,cosh_mass)
    observable_name = 'cosh mass state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(m_eff,m_eff_err))
    
def process_mass_log_jack_state_2(beta,N,SU,order,N_order,N_measure,N_thermal, accel = False):
    observable_name = 'state 2'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    measurement = np.load(file_name)
    stats = Stats.Stats(measurement)
    def log_mass(x):
    
        N_2 = int(N/2)
        X = np.concatenate((x,[x[0]]))
        cor = 1/2 * (X[:N_2+1] + X[N_2:][::-1])
        x_1 = np.roll(cor, -1)
        x_2 = np.roll(cor, 1)

        return - np.log(x_1/cor)
    m_eff,m_eff_err      =  stats.jacknife_arbitrary_function_2D_2_2(N,100,log_mass)
    observable_name = 'log mass state 2'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(m_eff,m_eff_err))

def process_Action_sq(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Action Sq'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    value, error,_,_ = Stats.Stats(data).estimate()
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(value,error))

def process_CV(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'Action Sq'

    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    action_sq,err_action_sq = np.load(file_name)
    observable_name = 'Action'

    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    action,err_action = np.load(file_name)
    print(err_action_sq)
    factor = 2*SU*N**2
    action = action * factor
    err_action = err_action * factor
    final_act_sq = action**2/factor
    final_err_act_sq = 2*action*err_action/factor
    print(final_err_act_sq)
    C_V = action_sq - final_act_sq
    err_C_V = np.sqrt(err_action_sq**2 + final_err_act_sq**2)
    observable_name = 'C_V'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(C_V,err_C_V))
def process_CV_jack(beta,N,SU,order,N_order,N_measure,N_thermal,accel = False):
    observable_name = 'CV_data'
    if accel == False:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    data = np.load(file_name)
    value,error = Stats.Stats(data).jacknife_C_V(100)
    factor = 2*SU*N**2

    value = value / factor
    error = error / factor
    print(value)

    observable_name = 'C_V'
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  +  " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  +  str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"

    np.save(file_name,(value,error))
    
