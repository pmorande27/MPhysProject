import unittest
import numpy as np
import matrix_routines as Mat
import scipy as sci
class TestMatrixRoutines(unittest.TestCase):
    def setUp(self) -> None:
        self.N = 10
        self.SU = 3
        self.Identity = np.zeros((self.N,self.N,self.SU,self.SU))
        for i in range(self.N):
            for j in range(self.N):
                self.Identity[i,j] = np.eye(self.SU)
        self.A = np.random.uniform(size = (self.N,self.N,self.SU,self.SU))
        
        self.C = np.random.uniform(size = (self.N,self.N,self.SU,self.SU))+ 1j*np.random.uniform(size = (self.N,self.N,self.SU,self.SU))
        self.two =  np.random.uniform(size = (self.N,self.N,2,2))+ 1j*np.random.uniform(size = (self.N,self.N,2,2))
    def test_multiply_matrices_zero_matrix(self):
        zero_matrix = np.zeros((self.N, self.N, self.SU, self.SU))
        expected_output = np.zeros((self.N, self.N, self.SU, self.SU))
        np.testing.assert_array_equal(Mat.multiply_matrices(self.A, zero_matrix), expected_output)

    def test_multiply_matrices_identity_matrix(self):
        expected_output = self.A
        np.testing.assert_array_equal(Mat.multiply_matrices(self.A, self.Identity), expected_output)

    def test_multiply_matrices_random_input(self):
        B = np.random.uniform(size=(self.N, self.N, self.SU, self.SU))
        # Compute the expected output using a known correct method
        expected_output = np.matmul(self.A, B)
        np.testing.assert_array_almost_equal(Mat.multiply_matrices(self.A, B), expected_output)
    def test_matrix_multiplication_with_identity(self):
        self.assertTrue(np.all(self.A ==Mat.multiply_matrices(self.A,self.Identity) ))
        self.assertTrue(np.all(Mat.multiply_matrices(self.A,self.Identity) ==Mat.multiply_matrices(self.A,self.Identity) ))
    
    def test_determinant(self):
        ones = np.ones((self.N,self.N))
        self.assertTrue(np.all(ones==Mat.determinant(self.Identity)))
    def test_determinant_zero_matrix(self):
        zero_matrix = np.array([[np.zeros((self.N, self.N))]])
        expected_output = [[0]]
        self.assertEqual(Mat.determinant(zero_matrix), expected_output)

    def test_determinant_diagonal_matrix(self):
        diagonal_matrix = np.array([[np.diag([1, 2, 3])]])
        expected_output = [[np.prod(np.diag(diagonal_matrix[0][0]))]]
        self.assertEqual(Mat.determinant(diagonal_matrix), expected_output)

    def test_determinant_singular_matrix(self):
        singular_matrix = np.array([[np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])]])
        expected_output = [[0]]
        self.assertEqual(Mat.determinant(singular_matrix), expected_output)
    def test_hermitian(self):
        B = np.zeros((self.N,self.N,self.SU,self.SU),complex)
        for i in range(self.N):
            for j in range(self.N):
                B[i,j] = np.conjugate(np.transpose(self.C[i,j]))
        self.assertTrue(np.all(B==Mat.dagger(self.C)))
    def test_exponential(self):
        D = Mat.exponential(self.C,1,3,1)
        self.assertTrue(np.all(D==self.Identity))
        M = Mat.exponential(self.C,2,3,1)
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
        D = Mat.exponential(F,0,2,0)

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


    def test_cross_product_known_input(self):
        a = np.array([1, 0, 3, 4])
        b = np.array([5, 6, 7, 8])
        c = np.array([9, 10, 11, 12])
        expected_output = [-8, 0, 24, -16]
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, b, c), expected_output)
    def test_cross_product_zero_input(self):
        a = np.zeros(4)
        b = np.zeros(4)
        c = np.zeros(4)
        expected_output = np.zeros(4)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, b, c), expected_output)
    def test_cross_product_negative_input(self):
        a = np.array([-1, 0, -3, -4])
        b = np.array([-5, -6, -7, -8])
        c = np.array([-9, -10, -11, -12])
        expected_output = -Mat.cross_product_in_4_d(a,b,c)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(-a, -b, -c), expected_output)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(-a, -b, c), Mat.cross_product_in_4_d(a,b,c))
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(-a, b, -c), Mat.cross_product_in_4_d(a,b,c))
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, -b, -c), Mat.cross_product_in_4_d(a,b,c))
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(-a, b, c), -Mat.cross_product_in_4_d(a,b,c))
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, -b, c), -Mat.cross_product_in_4_d(a,b,c))
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, b, -c), -Mat.cross_product_in_4_d(a,b,c))
    def test_cross_product_mixed_zero_input(self):
        a = np.array([-1, 0, -3, -4])
        b =  np.array([-5, -6, -7, -8])
        c = np.zeros(4)
        expected_output = np.zeros(4)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, b, c), expected_output)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(c, b, a), expected_output)
        np.testing.assert_array_equal(Mat.cross_product_in_4_d(a, c, b), expected_output)



            



    
