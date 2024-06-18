from Chiral import Chiral
import numpy as np
from Stats import Stats
import matplotlib.pyplot as plt
import matrix_routines as Mat
import Exceptions
import Stats
def exponential_measurements(beta, N, SU, order, N_order, N_measure,N_thermal,N_taus, Hot_start = True,accel =False):
    rates = np.zeros(len(N_taus))
    averages,errors = np.zeros(len(N_taus)),np.zeros(len(N_taus))
    averages_two,errors_two = np.zeros(len(N_taus)),np.zeros(len(N_taus))
    for i,N_tau in enumerate(N_taus):
        
        model = Chiral(N,beta,0,0,1,1/N_tau,N_tau,SU,order=order,order_N=N_order,Hot_start=True,accel=accel)
        model.Calibration_Runs(N_measure,N_thermal)
        rates[i] = model.accepted/model.tries
        print(rates[i])
        exp_H = np.exp(-model.DHc)
        averages[i],errors[i],_,_ = Stats.Stats(exp_H).estimate()
        averages_two[i],errors_two[i],_,_ = Stats.Stats(model.DHc).estimate()
        print(averages[i],errors[i])
    observable_name = "Exponential H"
    if accel == False:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    else:
        file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
    np.save(file_name,(N_taus,averages,errors,rates,averages_two,errors_two))
def calibration(beta, N, SU, order, N_order, N_tau_guess = 2, Hot_start = True,accel =False):
    N_tau = N_tau_guess
    print('Calibration with beta = ' + str(beta) + " N = " +str(N)+ " SU = " + str(SU) )
    up = 0.95
    low = 0.75
    max_count = 10
    results = [0 for i in range(max_count)]
    for i in range(max_count):
        epsilon = 1/N_tau
        calibration_runs = 10**3
        lat = Chiral(N, beta, 0,0,1,epsilon, N_tau, SU, order=order, order_N = N_order, Hot_start=Hot_start,accel=accel)
        lat.Calibration_Runs(calibration_runs, 1000)
        rate = lat.accepted/lat.tries
        d_rate = 0.75-rate
        results[i] = (rate-up,N_tau)
        print(rate,N_tau)
        

        new_N = int(np.rint(N_tau*(1+d_rate)))
        if rate <=up and rate >= low:
            if accel == False:
                file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order)
            else:
                file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order) + " Accel"
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
    if accel == False:
        file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order)
    else:
        file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order) + " Accel"
    np.save(file_name, [N_tau,1/N_tau])
    return N_tau
def load_calibration(beta, N, SU, order, N_order,accel = False):
    if accel == False:
        file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order)+'.npy'
    else:
        file_name = "ChiralParams/Chiral Calibration parameters beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+ " Order = " + str(order) + " N Order = " + str(N_order) + " Accel.npy"
    values = np.load(    file_name)
    return int(values[0]),values[1]
def lookup(d_rate,N_tau,results):
    for (x,y) in results:
        if abs(x) == d_rate and y == N_tau:
            return x
        
def measure(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name, Hot_start = True,accel =False, mass = 0.1):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            if accel == False:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            else:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order, Hot_start=Hot_start,accel=accel, mass = mass)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    #print(Stats(vals).estimate())
    np.save(file_name,results)
def measure_mass_dep(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name, mass,Hot_start = True,accel =True ):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Mass = "+str(mass)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order, Hot_start=Hot_start,accel=accel, mass = mass)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    #print(Stats(vals).estimate())
    np.save(file_name,results)
def measure_two_obs(beta, N, SU, order, N_order, N_measure,N_thermal, observable_1,observable_2, observable_name, Hot_start = True,accel =False, mass = 0.1):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            if accel == False:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            else:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order, Hot_start=Hot_start,accel=accel, mass = mass)
            results,rate = model.generate_measurements_two_obs(observable_1,observable_2)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    #print(Stats(vals).estimate())
    np.save(file_name,results)

def measure_func_2D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name,Hot_start = True,accel =False):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            if accel == False:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            else:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order,Hot_start=Hot_start,accel=accel)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1).swapaxes(1,2)

    np.save(file_name,vals)
def measure_func_1D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name, Hot_start = True,accel =False):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            if accel == False:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            else:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order,Hot_start=Hot_start,accel=accel)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1)
    np.save(file_name,vals)
def measure_func_4D(beta, N, SU, order, N_order, N_measure,N_thermal, observable, observable_name, Hot_start = True,accel =False):
    count = 0
    while True:
        try:
            if count == 10:    
                count = 0
                print('Recalibration')
                calibration(beta,N,SU,order,N_order,N_tau,Hot_start=Hot_start,accel=accel)
            if accel == False:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
            else:
                file_name = "ChiralResults/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+" Accel.npy"
            N_tau, epsilon = load_calibration(beta,N,SU,order, N_order,accel=accel)
            model = Chiral(N,beta,N_measure,N_thermal,1,epsilon,N_tau,SU,1,order=order, order_N=N_order,Hot_start=Hot_start,accel=accel)
            results,rate = model.generate_measurements(observable)
        except (Exceptions.ChiralExceptions):
            count+= 1
            continue
        break
    results = np.array(results)
    vals = results.swapaxes(0,1).swapaxes(1,2).swapaxes(2,3).swapaxes(3,4)
    print(vals.shape)

    np.save(file_name,vals)

