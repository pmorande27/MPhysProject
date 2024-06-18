#!/usr/bin/env python3
from argparse import ArgumentParser
import run

def main():
    #betas = [0.18,0.225,0.25,0.27,0.27,0.27,0.29,0.3,0.31,0.315]
    #Ns = [18,24,30,36,42,48,82,90,120,120]
    #betas = [0.3,0.31]
    #Ns = [90,100]

    ## SU 2
    #betas = [0.1,0.15,0.16,0.17]
    #Ns = [40 for b in betas]
    #betas2 = [0.18,0.19,0.21,0.22]
    #Ns2 = [64 for b in betas2]
    #betas3 = [0.23,0.24,0.25,0.26]
    #Ns3 = [96 for b in betas3]
    #betas4 = [0.27,0.28,0.29]
    #Ns4 = [160 for b in betas4]
    #betas5 = [0.3,0.31,0.32]
    #Ns5 = [224 for b in betas5]
    #betas = betas + betas2 + betas3 + betas4 + betas5
    #Ns = Ns + Ns2 + Ns3 + Ns4 + Ns5


    ## SU 6
    betas = [0.25,0.26,0.27,0.28,0.285,0.29,0.29,0.3,0.3,0.3,0.31,0.31,0.32,0.32]
    Ns = [18,24,24,30,30,30,36,36,42,48,54,60,76,82]


    #betas = [0.281,0.282,0.283,0.284,0.285,0.286,0.287,0.288,0.289,0.290,0.291,0.292,0.293,0.294,0.295,0.296,0.297,0.298,0.299]
    #Ns = [224 for b in betas]
    #betas2 = [0.291,0.292,0.293,0.294,0.295,0.296,0.297,0.298,0.299]
    #Ns2 = [90 for b in betas2]
    
    
    #betas3 = [0.301,0.302,0.303,0.304,0.305,0.306,0.307,0.308,0.309]
    #Ns3 = [100 for b in betas3]
    #betas4 = [0.310,0.311,0.312,0.313,0.314,0.315]
    #Ns4 = [110 for b in betas4]
    print(len(Ns))

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('array_pos', type=int, help='position in the array')
    args = parser.parse_args()
    i = args.array_pos
    N_measure = 10**5
    N_thermal = 10**4
    order = 10
    N_order = 10
    SU = 6
    accel = True
    beta = betas[i]
    N = Ns[i]
    #mass = mass_range[i]
    print("Initialising process with SU = %d, beta = %.3f, N = %d, accel %d"%(SU,beta,N,accel))
    #Chiral_run.calibration(beta,N,SU,order,N_order,10,Hot_start=True,accel=accel)
    #run.run_corr_length(beta,N,SU,order,N_order,N_measure,N_thermal,accel)
    run.run_iats(beta,N,SU,order,N_order,N_measure,N_thermal,accel)
    #run.run_state_two(beta,N,SU,order,N_order,N_measure,N_thermal,accel)
    #run.run_acceptance(beta,N,SU,order,N_order,N_measure,N_thermal,mass,accel)
    #run.run_cv(beta,N,SU,order,N_order,N_measure,N_thermal,accel)
    #run.run_action(beta,N,SU,order,N_order,N_measure,N_thermal,accel)
    print("End of process with SU = %d, beta = %.2f, N = %d, accel =  %s"%(SU,beta,N,accel))


if __name__ == '__main__':
    main()