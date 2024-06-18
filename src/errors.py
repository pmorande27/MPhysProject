import numpy as np
def errors(eg,ratio,erg,errratio):
    print(eg/ratio)
    err_w = np.sqrt((erg/ratio)**2+(eg*errratio/ratio**2)**2)
    print(err_w)
    print("e_w = %.2f \pm %.2f"%(eg/ratio,err_w))
errors( 3.027,0.986,0.014,0.005)