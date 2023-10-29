import numpy as np
import cupy as cp
import cProfile
import time
import Chiral
def main():
    a =Chiral.Chiral(16,0.1,10,100,1,1,2,2,1,5,5)
cProfile.run('main()')