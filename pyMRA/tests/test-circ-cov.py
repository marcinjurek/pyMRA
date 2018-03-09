import numpy as np
import pdb
import sys

sys.path.append('..')

import MRATools as mt
import matplotlib.pyplot as plt
from MRATools import dist

from sklearn.gaussian_process.kernels import Matern as skMatern



"""
This script tests an implementation of the ExpCovFunCric
function. It displays the covariance matrix and lets the 
user examine its pattern
"""


if __name__=='__main__':


    N = 101
    l=0.1
    sig=0.1
    
    locs = np.linspace(0, 1, N)[:-1].reshape((N-1, 1))
    pdb.set_trace()
    
    #S = mt.ExpCovFun(locs, circular=True, l=l)
    #mt.dispMat(S, title="Exponential")

    S = mt.Matern32(locs, circular=True, l=l)
    mt.dispMat(S, title="Matern 3/2")
    np.linalg.cholesky(S)
    
    S = mt.Matern52(locs, circular=True, l=l)
    mt.dispMat(S, title="Matern 5/2")
    np.linalg.cholesky(S)
    
    S = mt.GaussianCovFun(locs, circular=True, l=l)
    #mt.dispMat(S, title="Gaussian")
    #np.linalg.cholesky(S)
    
