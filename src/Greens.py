import numpy as np
import Stats
def Greens_mom(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values,errors  = np.load(file_name)
    Greens_zero_mom = np.zeros((N+1))
    error_zero = np.zeros((N+1))
    for i in range(N):
        Greens_zero_mom[i] = np.sum(values[i])
        error_zero[i] = np.sqrt(np.sum(errors[i]**2))
    Greens_zero_mom[N] = Greens_zero_mom[0]
    error_zero[N] = error_zero[0]
    observable_name = 'Greens 0 Mom'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    np.save(file_name,(Greens_zero_mom/Greens_zero_mom[0],error_zero/Greens_zero_mom[0]))
def susceptibility_from_complete_greens_function(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 2'
    file_name = "ChiralResults/Processed"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values,errors = np.load(file_name)
    sus = 0
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    #result += (i**2+j**2)*values[i][j]
                    sus += values[i][j][m][n]

    return sus/N**2
    
def second_moment_correletion_length_two(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values, error = np.load(file_name)
    result = 0
    sus = 0
    for i in range(N):
        for j in range(N):
            result += (i**2+j**2)*values[i][j]
            sus += values[i][j]

    return (result/(4*sus))**0.5,sus
def second_moment_correletion_length(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values, error = np.load(file_name)
    Gf = np.fft.fft2(values)
    p = 2*np.pi/N
    """Gf1 =np.sum(np.array([ [values[i,j] * np.exp(1j*(p)*j) for j in range(N)]for i in range(N)]))
    #print(np.sum(Gf1))
    Gf0 = np.sum(np.array([[values[i,j] for j in range(N)] for i in range(N)]))
    """
    print(Gf[0,1]-Gf[0,1].real)
    return (1/(4*np.sin(np.pi/N)*np.sin(np.pi/N))*(Gf[0,0]/Gf[0,1]-1))**0.5,Gf[0,0]

def second_moment_correletion_length_three(beta,N,SU,order,N_order,N_measure,N_thermal):
    observable_name = 'Greens 2'
    file_name = "ChiralResults/Processed/"+observable_name+"/"+observable_name+" beta = " + str(beta) + " N = " + str(N)  + " SU = " + str(SU)+" Order = "  + str(order)+" N Order = "  + str(N_order)+" N measurements = "  + str(N_measure)+" N Thermal = "  + str(N_thermal)+'.npy'
    values, error = np.load(file_name)
    result = 0
    sus = 0
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    result += ((i-m)**2+(j-n)**2)*values[i][j][m][n]
                    sus += values[i][j][m][n]
    result = result/N**2
    sus = sus /N**2
    result = result/(4*sus)
    return (result**0.5,sus)
