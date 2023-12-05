import numpy as np
import Stats

def process_Greens(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros((N,N))
    errors = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            values[i][j],errors[i][j],_ = Stats.Stats(data[i][j]).estimate()
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'

    np.save(file_name,(values,errors)) 
    
def process_ww_corr(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'ww corr'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'

    np.save(file_name,(values,values_err))
def process_G_0_mom(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 0 mom'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'

    np.save(file_name,(values,values_err))
def process_G_diags(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens diags'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros(N+1)
    values_err = np.zeros(N+1)
    for i in range(N):

        values[i],values_err[i],_= Stats.Stats(data[i]).estimate()
    values[N],values_err[N] = values[0],values_err[0]
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'

    np.save(file_name,(values,values_err))
def process_Greens_2(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 2'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    values = np.zeros((N,N,N,N))
    errors = np.zeros((N,N,N,N))
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    values[i][j][m][n],errors[i][j][m][n],_ = Stats.Stats(data[i][j][m][n]).estimate()
        
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    np.save(file_name,(values,errors))



def process_Action(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Action'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    value, error,_ = Stats.Stats(data).estimate()
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    np.save(file_name,(value,error))

def process_Generic_observable(beta,N,SU,order,N_order,N_measure,N_thermal, observable_name):
    
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    data = np.load(file_name)
    value, error,_ = Stats.Stats(data).estimate()
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    np.save(file_name,(value,error))

def reprocess_Greens(beta,N,SU,order,N_order,N_measure,N_thermal, N_thermal_new):
    observable_name = 'Greens'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    
    data = np.load(file_name)
    N_thermal_total = N_thermal + N_thermal_new
    values = np.zeros((N,N))
    errors = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            values[i][j],errors[i][j],_ = Stats.Stats(np.array(data[i][j])[N_thermal_new:]).estimate()
    N_measure = N_measure - N_thermal_new
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal_total)+'.npy'
    np.save(file_name,(values,errors)) 
def reprocess_ww(beta,N,SU,order,N_order,N_measure,N_thermal, N_thermal_new):
    observable_name = 'ww corr'
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    
    data = np.load(file_name)
    N_thermal_total = N_thermal + N_thermal_new
    values = np.zeros((N+1))
    errors = np.zeros((N+1))
    for i in range(N):
        values[i],errors[i],_ = Stats.Stats(np.array(data[i])[N_thermal_new:]).estimate()
    values[N],errors[N] = values[0],errors[0]
    N_measure = N_measure - N_thermal_new
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal_total)+'.npy'
    np.save(file_name,(values,errors)) 

    
