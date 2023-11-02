"""
Module used to perform all matrix operations
on the lattice
"""
import numpy as np

def dagger(lattice):
    """
    Method used to calculate the hermitian transpose
    of all the matrices in the lattice
    """
    return np.conjugate(np.einsum('ijkl->ijlk', lattice))

def multiply_matrices(lattice_a, lattice_b):
    """
    Method to calculate the matrix multiplication
    of all the elements (matrix_wise) of two lattices
    """

    return np.einsum('ijkl,ijlm->ijkm', lattice_a, lattice_b)

def exponential(lattice, order, su_parameter, order_n):
    """
    Method to calulcate the matrix exponential of all the matrices
    in the lattice to a given order.
    """
    if su_parameter == 2 and order == 0:

        generators = np.zeros((3, 2, 2), dtype=complex)

        generators[0][0, 1] = 1

        generators[0][1, 0] = 1

        generators[1][0, 1] = -1j

        generators[1][1, 0] = 1j

        generators[2][0, 0] = 1

        generators[2][1, 1] = -1

        p_is = np.einsum('ijkl,slk->ijs', lattice, generators).real/2

        result = exp_map(p_is)

        return result
    result = 0
    for n_order in range(order):
        #a = np.linalg.matrix_power(1j*U[:,:]/N,n)/np.math.factorial(n)

        result += np.linalg.matrix_power(1j * lattice[:, :]/order_n, n_order)\
                /np.math.factorial(n_order)

        #result+=  power(1j*U[:,:]/N,n,identity)/np.math.factorial(n)

    return np.linalg.matrix_power(result[:, :], order_n)

#   return power(result,N,identity)

def create_generators(su_parameter):
    """
    Method used to create the generators of SU(N)
    """
    generators = np.zeros((su_parameter**2-1, su_parameter, su_parameter), dtype=complex)

    if su_parameter == 3:

        generators[0][0, 1] = 1

        generators[0][1, 0] = 1

        generators[1][0, 1] = -1j

        generators[1][1, 0] = 1j

        generators[2][0, 0] = 1

        generators[2][1, 1] = -1

        generators[3][0, 2] = 1

        generators[3][2, 0] = 1

        generators[4][0, 2] = -1j

        generators[4][2, 0] = 1j

        generators[5][1, 2] = 1

        generators[5][2, 1] = 1

        generators[6][1, 2] = -1j

        generators[6][2, 1] = 1j

        generators[7][0, 0] = 1/np.sqrt(3)

        generators[7][1, 1] = 1/np.sqrt(3)

        generators[7][2, 2] = -2/np.sqrt(3)

    if su_parameter == 4:

        generators = create_generators_4()

    if su_parameter == 2:
        generators[0][0, 1] = 1

        generators[0][1, 0] = 1

        generators[1][0, 1] = -1j

        generators[1][1, 0] = 1j

        generators[2][0, 0] = 1

        generators[2][1, 1] = -1

    return generators

def create_generators_4():
    """
    Creates the generators for SU(4)
    """

    generators = np.zeros((4**2-1, 4, 4), dtype=complex)

    generators[0][0, 1] = 1

    generators[0][1, 0] = 1

    generators[1][0, 1] = -1j

    generators[1][1, 0] = 1j

    generators[2][0, 0] = 1

    generators[2][1, 1] = -1

    generators[3][0, 2] = 1

    generators[3][2, 0] = 1

    generators[4][0, 2] = -1j

    generators[4][2, 0] = 1j

    generators[5][1, 2] = 1

    generators[5][2, 1] = 1

    generators[6][1, 2] = -1j

    generators[6][2, 1] = 1j

    generators[7][0, 0] = 1/np.sqrt(3)

    generators[7][1, 1] = 1/np.sqrt(3)

    generators[7][2, 2] = -2/np.sqrt(3)

    generators[8][0, 3] = 1

    generators[8][3, 0] = 1

    generators[9][0, 3] = -1j

    generators[9][3, 0] = 1j

    generators[10][1, 3] = 1

    generators[10][3, 1] = 1

    generators[11][1, 3] = -1j

    generators[11][3, 1] = 1j

    generators[12][2, 3] = 1

    generators[12][3, 2] = 1

    generators[13][2, 3] = -1j

    generators[13][3, 2] = 1j

    generators[14][0, 0] = 1/np.sqrt(6)

    generators[14][1, 1] = 1/np.sqrt(6)

    generators[14][2, 2] = 1/np.sqrt(6)

    generators[14][3, 3] = -3/np.sqrt(6)

    return generators

def exp_map(alpha):
    """
    Method used to calculate the exponential map of the
    SU(2) group
    """
    n_length = alpha.shape[0]

    a_lattice = np.empty((n_length, n_length, 4))

    norm = np.sqrt(np.sum(alpha**2, axis=2)) # (N,N)

    # to do arithmetic with other (N,N,3) array need to broadcast to include axis 2

    alpha_norm = norm.reshape((n_length, n_length, 1))

    alpha_unit = np.divide(alpha, alpha_norm, out=np.zeros_like(alpha), where=alpha_norm != 0)

    a_lattice[:, :, 0] = np.cos(norm)

    a_lattice[:, :, 1:] = alpha_unit * np.sin(alpha_norm)

    mats = np.empty((n_length, n_length, 2, 2), dtype=complex)

    mats[:, :, 0, 0] = a_lattice[:, :, 0] + 1j * a_lattice[:, :, 3]

    mats[:, :, 0, 1] = a_lattice[:, :, 2]+ 1j * a_lattice[:, :, 1]

    mats[:, :, 1, 0] = -a_lattice[:, :, 2] + 1j * a_lattice[:, :, 1]

    mats[:, :, 1, 1] = a_lattice[:, :, 0] - 1j * a_lattice[:, :, 3]

    return mats

def determinant(lattice):
    """
    Methods used to calculate the determinant of all the matrices
    in the lattice
    """
    return np.linalg.det(lattice[:, :])


def reunitarisation(lattice, su_parameter):
    """
    Method used to reunitarise all the matrices in the
    lattice for SU(3)
    """
    if su_parameter == 3:

        n_length = len(lattice)

        for i in range(n_length):

            for j in range(n_length):

                lattice[i, j, 0] = lattice[i, j, 0]/np.linalg.norm(lattice[i, j, 0])

                lattice[i, j, 1] = lattice[i, j, 1] - \
                        np.dot(np.conjugate(lattice[i, j, 0]), lattice[i, j, 1])\
                        * lattice[i, j, 0]

                lattice[i, j, 1] = lattice[i, j, 1]/np.linalg.norm(lattice[i, j, 1])

                lattice[i, j, 2] = np.conjugate(np.cross((lattice[i, j, 0]), lattice[i, j, 1]))

    return lattice
