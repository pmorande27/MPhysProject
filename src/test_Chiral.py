import numpy as np
import unittest
from Chiral import Chiral

class TestChiral(unittest.TestCase):
    def setUp(self):
        self.chiral = Chiral(2, 1.0, 1, 0, 1, 0.1, 1, SU=3)

    def test_init(self):
        self.assertEqual(self.chiral.SU, 3)
        self.assertEqual(self.chiral.a, 1)
        self.assertEqual(self.chiral.order_N, 10)
        self.assertEqual(self.chiral.renorm_freq, 1000)
        self.assertEqual(self.chiral.epsilon, 0.1)
        self.assertEqual(self.chiral.N_tau, 1)
        self.assertEqual(self.chiral.N_sweeps, 1)
        self.assertEqual(self.chiral.c, 2)
        self.assertEqual(self.chiral.N, 2)
        self.assertEqual(self.chiral.N_measurement, 1)
        self.assertEqual(self.chiral.N_thermal, 0)
        self.assertEqual(self.chiral.beta, 1.0)

    def test_action(self):
        U = np.zeros((2, 2, 3, 3), dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    U[i, j, k, k] = 1
        beta = 1.0
        c = 2
        result = Chiral.action(U, beta, c)
        self.assertEqual(result, -24.0)

    def test_Hamiltonian(self):
        p = np.zeros((2, 2, 3, 3), dtype=complex)
        U = np.zeros((2, 2, 3, 3), dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    U[i, j, k, k] = 1
        beta = 1.0
        c = 2
        result = Chiral.Hamiltonian(p, U, beta, c)
        self.assertEqual(result, -24.0)

    
    def test_dot_p(self):
        U = np.zeros((2, 2, 3, 3), dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    U[i, j, k, k] = 1
        beta = 1.0
        SU = 3
        identity = U.copy()
        result = Chiral.dot_p(U, beta, SU, identity)
        np.testing.assert_array_almost_equal(result, np.zeros((2, 2, 3, 3), dtype=complex))

    def test_molecular_dynamics(self):
        p_0 = np.zeros((2, 2, 3, 3), dtype=complex)
        U_0 = np.zeros((2, 2, 3, 3), dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    U_0[i, j, k, k] = 1
        beta = 1.0
        epsilon = 0.1
        N_tau = 1
        SU = 3
        order = 10
        identity = U_0.copy()
        order_N = 10
        p, Unew = Chiral.molecular_dynamics(p_0, U_0, beta, epsilon, N_tau, SU, order, identity, order_N)
        np.testing.assert_array_almost_equal(p, np.zeros((2, 2, 3, 3), dtype=complex))
        np.testing.assert_array_almost_equal(Unew, U_0)

    def test_HMC(self):
        self.chiral.HMC()
        self.assertEqual(self.chiral.tries, 1)

    def test_thermalize(self):
        self.chiral.thermalize()
        self.assertEqual(self.chiral.tries, 0)

    def test_generate_measurements(self):
        self.chiral.U = np.random.normal(size =(2, 2, 3, 3))
        def observable(U):
            return 1.0
        results, rate = self.chiral.generate_measurements(observable)
        self.assertEqual(results, [1.0])
        self.assertEqual(rate, 1.0)

    def test_Calibration_Runs(self):
        self.chiral.Calibration_Runs(1, 1)
        self.assertEqual(self.chiral.tries, 1)

    def test_Measure_Sus(self):
        U = np.zeros((2, 2, 3, 3), dtype=complex)
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    U[i, j, k, k] = 1
        SU = 3
        result = Chiral.Measure_Sus(U, SU)
        self.assertEqual(result, 4.0)

   