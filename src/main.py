import numpy as np
import cProfile
import time
import Chiral
import Matrix_Routines as Mat
def main():
  U = np.random.random((10,10,3,3))
  A = Mat.reunitarisation(U,3)
  print(Mat.multiply_matrices(U,Mat.dagger(A)))
  print(np.linalg.det(Mat.reunitarisation(U,3)))
main()