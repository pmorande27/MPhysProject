from Chiral import Chiral
import numpy as np
from Stats import Stats
import matplotlib.pyplot as plt
import matrix_routines as Mat
import Exceptions

def calibration(beta, N, SU, order, N_order, N_tau_guess = 2, Hot_start = True):
    N_tau = N_tau_guess
    print('Calibration with beta = ' + str(beta) + " N = " +str(N)+ " SU = " + str(SU) )
    up = 0.95
    low = 0.75
    max_count = 10
    results = [0 for i in range(max_count)]
    for i in range(max_count):
        epsilon = 1/N_tau
        calibration_runs = 10**3
        lat = Chiral(N, beta, 0,0,1,epsilon, N_tau, SU, order=order, order_N = N_order, Hot_start=Hot_start)
        lat.Calibration_Runs(calibration_runs, 1000)
        rate = lat.accepted/lat.tries
        d_rate = 0.75-rate
        results[i] = (rate-up,N_tau)
        print(rate,N_tau)
        

        new_N = int(np.rint(N_tau*(1+d_rate)))
        if rate <=up and rate >= low:
            file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order)
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
    file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order)

    np.save(file_name, [N_tau,1/N_tau])
    return N_tau
def load_calibration(beta, N, SU, order, N_order):
    file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order) +'.npy'
    values = np.load(    file_name)
    return int(values[0]),values[1]
def lookup(d_rate,N_tau,results):
    for (x,y) in results:
        if abs(x) == d_rate and y == N_tau:
            return x
        
def measure(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau)
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    #print(Stats(vals).estimate())
    np.save(file_name,results)

def measure_func_2D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau)
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1).swapaxes(1,2)

    np.save(file_name,vals)
def measure_func_1D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name, Hot_start = True):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau)
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order,Hot_start=Hot_start)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1)
    np.save(file_name,vals)
def measure_func_4D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau)
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1).swapaxes(1,2).swapaxes(2,3).swapaxes(3,4)
    print(vals.shape)

    np.save(file_name,vals)

def measure_susceptibility(beta,N,SU):
    file_name = "ChiralResults/Susceptibility/Susceptibility 2 beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+'.npy'
    N_tau, epsilon = load_calibration(beta,N,SU)
    model = Chiral(N,beta,10**3,1000,1,epsilon,N_tau,2,1,order=10)
        
    results = model.generate_measurements(lambda U:np.einsum('ijkl,ablk->',U,Mat.dagger(U)).real/(2*N**2))
    np.save(file_name,results)


def estimate_2D_func(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name):
    file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    vals = np.load(file_name)
    results = np.zeros((N,N))
    err = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            a = Stats(vals[i,j]).estimate()
            results[i,j] = a[0]
            err[i,j] = a[1]
    return results
