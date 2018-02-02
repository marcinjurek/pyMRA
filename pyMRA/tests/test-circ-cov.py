import numpy as np
import pdb
import sys

sys.path.append('..')

import MRATools as mt
import matplotlib.pyplot as plt



"""
This script tests an implementation of the ExpCovFunCric
function. It displays the covariance matrix and lets the 
user examine its pattern
"""


if __name__=='__main__':


    N1 = 10; N2 = 15
        
    locs1 = np.linspace(0, 1, N1)[:-1]#.reshape((N1-1, 1))
    locs2 = np.linspace(0, 1, N2)[:-1]#.reshape((N2-1, 1))
    S = mt.ExpCovFun(locs1, circular=True)
    mt.dispMat(S)




    
    locs = np.linspace(start=0, stop=1, num=11)[:-1]
    #locs = np.array([0, 0.25, 0.5, 0.75])
    
    xv, yv = np.meshgrid(locs, locs)
    m = np.minimum(xv, yv)
    M = np.maximum(xv, yv)
    dist = np.minimum(M - m, m + 1-M)
    dist2 = np.roll(dist, 3, axis=1)
   
    S = mt.ExpCovFunCirc(locs)
    mt.dispMat(S)
