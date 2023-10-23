import numpy as np

def dagger(U):
    return np.conjugate(np.einsum('ijkl->ijlk',U))
def multiply_matrices(A,B):
     
    return np.einsum('ijkl,ijlm->ijkm',A,B)

def power(U,n,identity):
    N = len(U)
    SU = len(U[0][0])
    result = identity.copy()

    if n == 0:
        return result
    else:
        for i in range(n):
            result = np.einsum('ijkl,ijlm->ijkm',result,U)
        return result

def exponential(U, order,SU,identity,N):
    N = len(U)
    """if SU == 2:
        generators = np.zeros((3,2,2),dtype=complex)

        generators[0][0,1] =1
        generators[0][1,0] = 1

        generators[1][0,1] = -1j
        generators[1][1,0] = 1j

        generators[2][0,0] = 1
        generators[2][1,1] =-1
        
        p_is  = np.einsum('ijkl,slk->ijs',U,generators).real/2
        result = exp_map(p_is)
        return result"""
    result = 0
    for n in range(order):
        #a = np.linalg.matrix_power(1j*U[:,:]/N,n)/np.math.factorial(n)
        result+=  np.linalg.matrix_power(1j*U[:,:]/N,n)/np.math.factorial(n)
        #result+=  power(1j*U[:,:]/N,n,identity)/np.math.factorial(n)

    
    return np.linalg.matrix_power(result[:,:],N)
#   return power(result,N,identity)

def create_generators(SU):
    generators = np.zeros((SU**2-1,SU,SU),dtype=complex)

    if SU == 3:
        generators[0][0,1] = 1
        generators[0][1,0] = 1
        generators[1][0,1] = -1j
        generators[1][1,0] = 1j
        generators[2][0,0]  = 1
        generators[2][1,1] = -1
        generators[3][0,2] = 1
        generators[3][2,0] = 1
        generators[4][0,2] = -1j
        generators[4][2,0] = 1j
        generators[5][1,2] = 1
        generators[5][2,1] = 1
        generators[6][1,2] = -1j
        generators[6][2,1] = 1j
        generators[7][0,0] = 1/np.sqrt(3)
        generators[7][1,1] = 1/np.sqrt(3)
        generators[7][2,2] = -2/np.sqrt(3)
    if SU == 2:
        generators[0][0,1] =1
        generators[0][1,0] = 1

        generators[1][0,1] = -1j
        generators[1][1,0] = 1j

        generators[2][0,0] = 1
        generators[2][1,1] =-1
    return generators

def exp_map(alpha):
    N = alpha.shape[0]
    a = np.empty((N,N,4))
    norm = np.sqrt(np.sum(alpha**2, axis=2)) # (N,N)
    # to do arithmetic with other (N,N,3) array need to broadcast to include axis 2
    alpha_norm = norm.reshape((N,N,1))
    alpha_unit = np.divide(alpha, alpha_norm, out=np.zeros_like(alpha), where=alpha_norm!=0) # avoids division by zero. When norm is zero, i.e alpha is zero, alpha_unit is set to zero too 
    a[:,:,0] = np.cos(norm)
    a[:,:,1:] = alpha_unit * np.sin(alpha_norm)
    mats = np.empty((N,N,2,2), dtype=complex)
    mats[:,:,0,0] = a[:,:,0]+1j*a[:,:,3]
    mats[:,:,0,1] = a[:,:,2]+1j*a[:,:,1]
    mats[:,:,1,0] = -a[:,:,2]+1j*a[:,:,1]
    mats[:,:,1,1] = a[:,:,0]-1j*a[:,:,3]

    
    return mats
   
