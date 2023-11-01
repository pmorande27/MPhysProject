import unittest
import numpy as np
import Matrix_Routines as Mat
import scipy as sci
class TestSolarSystem(unittest.TestCase):
    def setUp(self) -> None:
        self.N = 2
        self.SU = 3
        self.Identity = np.zeros((self.N,self.N,self.SU,self.SU))
        for i in range(self.N):
            for j in range(self.N):
                self.Identity[i,j] = np.eye(self.SU)
        self.A = np.random.uniform(size = (self.N,self.N,self.SU,self.SU))
        
        self.C = np.random.uniform(size = (self.N,self.N,self.SU,self.SU))+ 1j*np.random.uniform(size = (self.N,self.N,self.SU,self.SU))
        self.two =  np.random.uniform(size = (self.N,self.N,2,2))+ 1j*np.random.uniform(size = (self.N,self.N,2,2))
    def test_matrix_multiplication_with_identity(self):
        self.assertTrue(np.all(self.A ==Mat.multiply_matrices(self.A,self.Identity) ))
        self.assertTrue(np.all(Mat.multiply_matrices(self.A,self.Identity) ==Mat.multiply_matrices(self.A,self.Identity) ))
    
    def test_matrix_power(self):
        self.assertTrue(np.all(Mat.multiply_matrices(self.A,self.A)==Mat.power(self.A,2,self.Identity)))
        self.assertTrue(np.all(self.Identity==Mat.power(self.A,0,self.Identity)))
    def test_determinant(self):
        ones = np.ones((self.N,self.N))
        self.assertTrue(np.all(ones==Mat.determinant(self.Identity)))
    
    def test_hermitian(self):
        B = np.zeros((self.N,self.N,self.SU,self.SU),complex)
        for i in range(self.N):
            for j in range(self.N):
                B[i,j] = np.conjugate(np.transpose(self.C[i,j]))
        self.assertTrue(np.all(B==Mat.dagger(self.C)))
    def test_exponential(self):
        D = Mat.exponential(self.C,1,3,self.Identity,1)
        self.assertTrue(np.all(D==self.Identity))
        M = Mat.exponential(self.C,2,3,self.Identity,1)
        self.assertTrue(np.all(M==(self.Identity+1j*self.C)))
    
    def test_reunitaristation(self):
        self.assertTrue(np.all(self.Identity==Mat.reunitarisation(self.Identity,3)))
        np.testing.assert_almost_equal(np.ones((self.N,self.N)),Mat.determinant(Mat.reunitarisation(self.C,3)),decimal=15)
        np.testing.assert_almost_equal(self.Identity,Mat.multiply_matrices(self.C,Mat.dagger(Mat.reunitarisation(self.C,3))),decimal=15)

    def test_exponential_two(self):
        
        generators = np.zeros((3,2,2),complex)
        generators[0][0,1] =1
        generators[0][1,0] = 1

        generators[1][0,1] = -1j
        generators[1][1,0] = 1j

        generators[2][0,0] = 1
        generators[2][1,1] =-1
        a = np.random.random((16,16,3))
        F =np.einsum('ijk,klm->ijlm',a,generators)
        D = Mat.exponential(F,0,2,self.Identity,0)
        E = np.zeros((self.N,self.N,2,2),complex)
        for i in range(self.N):
            for J in range(self.N):
                E = sci.linalg.expm(1j*F)
        np.testing.assert_almost_equal(D,E,15)

    def test_generators(self):
        SUs = [2,3,4]
        c = [2.,2.,2.]
        for i,SU in enumerate(SUs):
            generators = Mat.create_generators(SU)
            trace = np.zeros(len(generators))
            np.testing.assert_equal(trace,np.einsum('ikk',generators))
            np.testing.assert_equal(generators,np.conjugate(np.einsum('ijk->ikj',generators)))
            np.testing.assert_almost_equal(np.einsum('ijk,lkj',generators,generators),c[i]*(np.eye(SU**2-1)+np.zeros((SU**2-1,SU**2-1),complex)))



            



    
